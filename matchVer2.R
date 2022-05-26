
library(Seurat)
library(HandyPack)
library(stringr)

rm(list=ls())

source('CopyWare.R')
source('IntegrationWare.R')
## ####################################################
## ####################################################

## ####################################################
## Subset the raw objects:
oldCombined = getOld('combined')
oldCombined = tagCells(oldCombined)

## Raw objects:
tags = c('2dpf','3dpf','adult19','adult20')

obj = list()
copyQC = list()

for(tag in tags)
{
    f = readRDS(paste0('unintegrated/',tag,'.rds'))
    f = tagCells(f)
    obj[[tag]] = f

    idx = f$tag %in% oldCombined$tag
    f = f[,idx]
    copyQC[[tag]] = f

    fileOut = paste0('unintegrated/',tag,'_copyQC.rds')
    saveRDS(f,fileOut)
}

## ####################################################
## Make the per-stage objects:

## Make the integrated larval object:
day2 = NormalizeData(copyQC[['2dpf']])
day3 = NormalizeData(copyQC[['3dpf']])

day2 = FindVariableFeatures(day2)
day3 = FindVariableFeatures(day3)

anchors = FindIntegrationAnchors(list(day2,day3))
larval = IntegrateData(anchors)

larval = runBasicAnalyses(larval)  

saveDir = nameAndMakeDir('intermediateCopyQC')
saveRDS(larval,
        'intermediateCopyQC/larval_integrated.rds')

## Make the adult objects:
adult19 = copyQC[['adult19']]
adult20 = copyQC[['adult20']]

## 20 only:
adult20Only = NormalizeData(adult20)
adult20Only = runBasicAnalyses(adult20Only)
saveRDS(adult20Only,'intermediateCopyQC/adult20Only.rds')

## Merged:
adult = merge(adult19,y=adult20,
              project='audlt')
adult = NormalizeData(adult)
adult = runBasicAnalyses(adult)
saveRDS(adult,'intermediateCopyQC/adult.rds')



## ####################################################
## Copy CellType into objects:

## Larval:
oldLarval = getOld('larval')
oldLarval = tagCells(oldLarval)
dictionary = oldLarval$CellType
names(dictionary) = oldLarval$tag
larval$CellType = dictionary[larval$tag]
larval$CellType[is.na(larval$CellType)]  = '*'

saveRDS(larval,
        'intermediateCopyQC/larval_CellType.rds')

## Adult:
oldAdult = getOld('adult')
oldAdult = tagCells(oldAdult)
dictionary = oldAdult$shortName
names(dictionary) = oldAdult$tag

## 20 only:
adult20Only$CellType = dictionary[adult20Only$tag]
adult20Only$CellType[is.na(adult20Only$CellType)]  = '*'
saveRDS(adult20Only,
            'intermediateCopyQC/adult20Only_CellType.rds')

## All:
adult$CellType = dictionary[adult$tag]
adult$CellType[is.na(adult$CellType)]  = '*'
saveRDS(adult,
        'intermediateCopyQC/adult_CellType.rds')

## ####################################################
## Make the combined object:

## This is still going to vary by keepCells:
for(whichKeep in 1:2)
{
    keepCellTypesDF = Read.Table('keepCellTypes.txt')
    if(whichKeep == 1)
    {
        idx = keepCellTypesDF$Keep == 'Y'
    } else {
        idx = keepCellTypesDF$Keep2 == 'Y'
    }
    keepCellTypes = keepCellTypesDF$CellType[idx]

    idx = larval$CellType %in% keepCellTypes
    larval = larval[,idx]

    ## Combined adult20Only:
    fileOut = paste0('intermediateCopyQC/combined20Only_keepCells',
                     whichKeep,
                     '.rds')
    anchors = FindIntegrationAnchors(list(adult20Only,larval))
    combined20Only = IntegratData(anchors)
    combined20Only = runBasicAnalyses(combined20Only)
    saveRDS(combined20Only,fileOut)
    
    ## All
    fileOut = paste0('intermediateCopyQC/combined_keepCells',
                     whichKeep,
                     '.rds')
    anchors = FindIntegrationAnchors(list(adult,larval))
    combined = IntegratData(anchors)
    combined = runBasicAnalyses(combined)
    saveRDS(combined,fileOut)
}

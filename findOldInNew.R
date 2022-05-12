
library(Seurat)
library(stringr)

rm(list=ls())

source('ZebraWare.R')
## ####################################################
## ####################################################
tagCells = function(obj)
{
    obj = stripCellNames(obj)
    obj$tag = paste(obj$stripped,
                    obj$orig.ident,
                    sep='_')

    return(obj)
}

## ####################################################
stripCellNames = function(obj)
{
    headers = c('adult19_adult19_',
                'adult19_',
                'adult20_adult20_',
                'adult20_',
                '2dpf_2dpf_',
                '2dpf_')
                
    n = colnames(obj)
    for(h in headers)
        n = str_replace(n,h,'')
 
    n = str_replace(n,'-.*','')

    obj$stripped = n
    return(obj)
}

## ####################################################
getOld = function(which)
{
    if(which == 'combined')
        fileName = '~/TiffanyZebraFish02/SeuratObject/combinedSubclusteredNamed.rds'

    if(which == 'larval')
        fileName = paste0('~/TiffanyZebraFish/larvalZebraFish/data/',
                          'SeuratObject.rds')

    if(which == 'adult')
        fileName = '~/TiffanyZebraFish02/SeuratObject/adult.rds'
    
    f = readRDS(fileName)
    if('orig_ident' %in% names(f@meta.data))
        f$orig.ident = f$orig_ident

    

    ## f = uniformOrig(f)

    return(f)
}

## ####################################################
uniformOrig = function(f)
{
    keys = c('adult19','adult20')
    for(k in keys)
    {
        idx = str_detect(f$orig.ident,k)
        f$orig.ident = k
    }

    return(f)
}

## ####################################################
orig = function(which)
{
    orig = c(old='orig.ident')

    return(orig[which])
}

## ####################################################
testDuplication = function(obj)
{
    return(sum(duplicated(obj$tag)) > 0)
}

## ####################################################
old = list()
old$larval = getOld('larval')
old$adult = getOld('adult')
old$combined = getOld('combined')

new = list()
new$larval = getSeuratObject('larval')
new$adult = getSeuratObject('adult')
new$combined = getSeuratObject('combined')

writeLines('old')
for(n in names(old))
{
    writeLines(n)
    writeLines(as.character(unique(old[[n]]$orig.ident)))
    writeLines('')
}

writeLines('new')
for(n in names(new))
{
    writeLines(n)
    writeLines(as.character(unique(new[[n]]$orig.ident)))
    writeLines('')
}

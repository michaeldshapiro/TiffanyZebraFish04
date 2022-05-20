
library(Seurat)
library(HandyPack)
library(stringr)

## ####################################################
## ####################################################
enforceOrig.ident = function(f)
{
    idx = names(f@meta.data) == 'orig_ident'
    names(f@meta.data)[idx] = 'orig.ident'

    return(f)
}

## ####################################################
uniformOrig.ident = function(f)
{
    dictionary = c(adult191217='adult19',
                   adult200108='adult20',
                   larval48='2dpf',
                   larval68='3dpf',
                   HEA614A1='adult20',
                   '48HPF'='2dpf',
                   '69HPF'='3dpf')

    for(n in names(dictionary))
        f$orig.ident = str_replace(f$orig.ident,n,dictionary[n])

    return(f)
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

    f = uniformOrig.ident(f)

    return(f)
}
 
## ####################################################
stripCellNames = function(obj)
{
    headers = c('adult19_',
                'adult20_',
                '2dpf_',
                '3dpf_')
                
    n = colnames(obj)
    for(h in headers)
        n = str_replace(n,h,'')
 
    n = str_replace(n,'-.*','')
    n = str_replace(n,'_.*','')

    obj$stripped = n
    return(obj)
}

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
copyCellTypeIntoNewLarval = function()
{
    larval = readRDS('intermediate/larval.rds')
    oldLarval = getOld('larval')

    larval = tagCells(larval)
    oldLarval = tagCells(oldLarval)

    dictionary = oldLarval$CellType
    names(dictionary) = oldLarval$tag

    larval$CellType = dictionary[larval$tag]
    larval$CellType[is.na(larval$CellType)]  = '*'
    
    saveRDS(larval,
            'intermediate/larval_CellType.rds')
}

## ####################################################
copyCellTypeIntoNewIntegratedLarval = function()
{
    larval = readRDS('intermediate/larval_integrated.rds')
    oldLarval = getOld('larval')

    larval = tagCells(larval)
    oldLarval = tagCells(oldLarval)

    dictionary = oldLarval$CellType
    names(dictionary) = oldLarval$tag

    larval$CellType = dictionary[larval$tag]
    larval$CellType[is.na(larval$CellType)]  = '*'
    
    saveRDS(larval,
            'intermediate/larval_integrated_CellType.rds')
}


## ####################################################
copyCellTypeIntoNewAdult = function()
{
    adult = readRDS('intermediate/adult.rds')
    oldAdult = getOld('adult')

    adult = tagCells(adult)
    oldAdult = tagCells(oldAdult)

    dictionary = oldAdult$shortName
    names(dictionary) = oldAdult$tag

    adult$CellType = dictionary[adult$tag]
    adult$CellType[is.na(adult$CellType)]  = '*'
    
    saveRDS(adult,
            'intermediate/adult_CellType.rds')
}


## ####################################################
copyCellTypeIntoNewAdult20 = function()
{
    adult = readRDS('unintegrated/adult20_qc.rds')
    oldAdult = getOld('adult')

    adult = tagCells(adult)
    oldAdult = tagCells(oldAdult)

    dictionary = oldAdult$shortName
    names(dictionary) = oldAdult$tag

    adult$CellType = dictionary[adult$tag]
    adult$CellType[is.na(adult$CellType)]  = '*'
    
    saveRDS(adult,
            'unintegrated/adult20_CellType.rds')
}


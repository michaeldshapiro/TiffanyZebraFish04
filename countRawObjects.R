
library(Seurat)
library(HandyPack)
library(ggplot2)
library(stringr)

rm(list=ls())
graphics.off()

## ####################################################
## ####################################################
getObject = function(version,whichOne)
{
    fileList = getFiles()
    
    if(version == 1)
    {
        f = readRDS(fileList[[1]])
        f = extractSubset(f,whichOne)
        if(is.null(f))
            return(f)
    } else {
        f = readRDS(fileList[[version]][whichOne])
    }

    if(! 'percent.mt' %in% names(f@meta.data))
    {
        b = PercentageFeatureSet(f, pattern = "^mt-")
        f$percent.mt = b$nCount_RNA
    }

    if('orig_ident' %in% names(f@meta.data))
        f$orig.ident = f$orig_ident

    return(f)
}

## ####################################################
extractSubset = function(f,whichOne)
{
    if(whichOne == 'adult19')
        return(NULL)
    
    tr = c('2dpf'='sc50hpf',
           '3dpf'='sc70hpf',
           adult20='snAdult',
           adult19='')

    fThis = f[,f$orig_ident==tr[whichOne]]

    return(fThis)    
}

## ####################################################
getFiles = function()
{
    files = list()

    ## version 1 is spurious:
    files[[1]] = paste0('~/TiffanyZebraFish/',
                        'SC18139_and_GSE_152906/data/provisional.rds')
    ## version 2:
    files[[2]] = c(adult19='~/TiffanyZebraFish02/Integration01/unintegrated/adult191217.rds',
                   adult20='~/TiffanyZebraFish02/Integration01/unintegrated/adult200108.rds',
                   '2dpf'='~/TiffanyZebraFish02/Integration01/unintegrated/larval48.rds',
                   '3dpf'='~/TiffanyZebraFish02/Integration01/unintegrated/larval68.rds')

    ## version 3:
    files[[3]] = c(adult19='~/TiffanyZebraFish03/Integration/unintegrated/adult19.rds',
                   adult20='~/TiffanyZebraFish03/Integration/unintegrated/adult20.rds',
                   '2dpf'='~/TiffanyZebraFish03/Integration/unintegrated/2dpf.rds',
                   '3dpf'='~/TiffanyZebraFish03/Integration/unintegrated/3dpf.rds')

    ## version 4:
    files[[4]] = c(adult19='~/TiffanyZebraFish04/unintegrated/adult19.rds',
                   adult20='~/TiffanyZebraFish04/unintegrated/adult20.rds',
                   '2dpf'='~/TiffanyZebraFish04/unintegrated/2dpf.rds',
                   '3dpf'='~/TiffanyZebraFish04/unintegrated/3dpf.rds')

    return(files)
}

## ####################################################
## ####################################################
idents = c('2dpf','3dpf','adult19','adult20')
versions = 1:4

vars = c('percent.mt','nFeature_RNA','nCount_RNA')

finger = 1

for(version in versions)
{
    for(ident in idents)
    {
        f = getObject(version,ident)
        if(is.null(f))
            next
        a = FetchData(f,vars)
        a$version = paste0('version_',version)
        a$stage = ident 

        if(finger == 1)
            df = a
        else
            df = rbind(df,a)

        finger = finger + 1
    }
}

for(var in vars)
{
    for(stage in idents)
    {
        idx = df$stage == stage

        title = paste(var,'in',stage)
        g = ggplot(df[idx,],aes_string(x=var)) +
            geom_histogram(bins=50) +
            ggtitle(title) +
            facet_wrap(~version,ncol=1)

        title = str_replace(title,' ','_')
        saveDir = nameAndMakeDir(c('counts','rawObjects'))
        fileName = paste0(saveDir,'/',title,'.jpg')
        dev.new()
        print(g)
        ggsave(plot=g,
               filename=fileName,
               height=8,width=14,units='in')
    }
}

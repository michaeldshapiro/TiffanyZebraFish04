
library(Seurat)
library(HandyPack)
library(ggplot2)
library(stringr)

rm(list=ls())
graphics.off()

## ####################################################
## ####################################################
getCombined = function(which)
{
    if(which == 1)
    {
        fileName =  paste0('~/TiffanyZebraFish/',
                           'SC18139_and_GSE_152906/data/provisional.rds')
        f = readRDS(fileName)
        f$orig.ident = f$orig_ident
        f$percent.mt = f$percent_mt
        return(f)
    }

    if(which == 2)
    {
        f = readRDS('~/TiffanyZebraFish02/SeuratObject/combinedTSNE.rds')
        return(f)
    }

    if(which == 4)
    {
        fileName = paste0('~/TiffanyZebraFish04/intermediate/',
                          'combined_integratedLarval_keepCells2_adult20.rds')
        f = readRDS(fileName)
        return(f)
    }
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
                   '69HPF'='3dpf',
                   sc50hpf='2dpf',
                   sc70hpf='3dpf',
                   snAdult='adult20')

    for(n in names(dictionary))
        f$orig.ident = str_replace(f$orig.ident,n,dictionary[n])

    return(f)
}

## ####################################################
cells = c()
cells[1] = 'shorterName'
cells[2] = 'cellType'
cells[4] = 'CellType'
combined = list()

## ####################################################
## We want versions 1,2,4:
versions = c(1,2,4)
for(i in versions)
    combined[[i]] = uniformOrig.ident(getCombined(i))

idents = c('2dpf','3dpf','adult19','adult20','day5')
M = matrix(0,nrow=3,ncol=5)
rownames(M) = paste0('version_',versions)
colnames(M) = idents


for(i in 1:3)
    for(j in 1:5)
        M[i,j] = sum(combined[[versions[i]]]$orig.ident == idents[j])
        
df = data.frame(M)
names(df) = idents

stop('look and listen')

saveDir = nameAndMakeDir('counts') 
Write.Table(df,'counts/orig.identCounts.txt')

## ####################################################
## Percent mt is only in two of these:
variables = c('percent.mt','nFeature_RNA','nCount_RNA')
for(i in versions[c(1,3)])
{
    ## larval:
    idx = combined[[i]]$orig.ident %in% idents[1:2]
    fLarval = combined[[i]][,idx]
    aLarval = FetchData(fLarval,variables)
    aLarval$stage = 'larval'
    
    ## adult:
    idx = combined[[i]]$orig.ident %in% idents[3:4]
    fAdult = combined[[i]][,idx]
    aAdult = FetchData(fAdult,variables)
    aAdult$stage = 'adult'

    a = rbind(aLarval,aAdult)
    a$version = paste0('version_',i)

    if(i == 1)
        cutoffsDF = a
    else
        cutoffsDF = rbind(cutoffsDF,a)
}

g = ggplot(cutoffsDF,aes(x=percent.mt)) +
    geom_histogram(bins=50) +
    facet_wrap(~ version + stage) +
    ggtitle('Percent Mt')

dev.new()
print(g)
GGSave(plot=g,
       filename='counts/combined_percentMT.jpg',
       width=12)

## ####################################################
## Percent mt is only in two of these:
variables = c('nFeature_RNA','nCount_RNA')
for(i in versions)
{
    ## larval:
    idx = combined[[i]]$orig.ident %in% idents[1:2]
    fLarval = combined[[i]][,idx]
    aLarval = FetchData(fLarval,variables)
    aLarval$stage = 'larval'
    
    ## adult:
    idx = combined[[i]]$orig.ident %in% idents[3:4]
    fAdult = combined[[i]][,idx]
    aAdult = FetchData(fAdult,variables)
    aAdult$stage = 'adult'

    a = rbind(aLarval,aAdult)
    a$version = paste0('version_',i)

    if(i == 1)
        cutoffsDF = a
    else
        cutoffsDF = rbind(cutoffsDF,a)
}

h = ggplot(cutoffsDF,aes(x=nCount_RNA)) +
    geom_histogram(bins=50) +
    facet_wrap(~ version + stage,ncol=2) +
    ggtitle('nCount_RNA')

dev.new()
print(h)
GGSave(plot=h,
       filename='counts/combined_nCount_RNA.jpg',
       width=12)


k = ggplot(cutoffsDF,aes(x=nFeature_RNA)) +
    geom_histogram(bins=50) +
    facet_wrap(~ version + stage,ncol=2) +
    ggtitle('nFeature_RNA')

dev.new()
print(k)
GGSave(plot=k,
       filename='counts/combined_nFeature_RNA.jpg',
       width=12)

## ####################################################
## Extract min:
theVersions = unique(cutoffsDF$version)
minLarvalCount = c()
minAdultCount = c()
minLarvalFeature = c()
minAdultFeature = c()

larval = cutoffsDF[cutoffsDF$stage=='larval',]
adult = cutoffsDF[cutoffsDF$stage=='adult',]

for(version in theVersions)
{
    minLarvalCount[version] = min(larval[larval$version==version,]$nCount_RNA)
    minLarvalFeature[version] = min(larval[larval$version==version,]$nFeature_RNA)

    minAdultCount[version] = min(adult[adult$version==version,]$nCount_RNA)
    minAdultFeature[version] = min(adult[adult$version==version,]$nFeature_RNA)    
}

minima = data.frame(version=theVersions,
                    minLarvalCount,
                    minLarvalFeature,
                    minAdultCount,
                    minAdultFeature)
Write.Table(minima,
            'counts/minFeaturesAndCounts.txt')


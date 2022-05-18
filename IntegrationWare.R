
library(Seurat)
library(HandyPack)
library(tictoc)
library(stringr)
library(ggplot2)
library(plotly)

## ####################################################
## ####################################################
##
## 04
##
## ####################################################
## ####################################################
makeSeuratObjectsFrom10X = function()
{
    dataDirs = c('2dpf'='data/2dpf',
                 '3dpf'='data/3dpf',
                 adult19='data/adult/191217',
                 adult20='data/adult/200108')


    outDir = nameAndMakeDir('unintegrated')

    for(n in names(dataDirs))
    {
        Tic(n)
        where = dataDirs[n]
        obj = Read10X(where)
        SeuratObject = CreateSeuratObject(counts=obj,
                                          project=n,
                                          min.cells=10)

        ## Add the mt data:
        b = PercentageFeatureSet(SeuratObject, pattern = "^mt-")
        SeuratObject$percent.mt = b$nCount_RNA
        fileName = paste0(outDir,
                          '/',
                          n,
                          '.rds')

        saveRDS(SeuratObject,fileName)
        toc()
    }
}

## ####################################################
makeQCFigures = function()
{
    figDir = nameAndMakeDir('integrationFigures')
    files = Sys.glob('unintegrated/*')
    idx = str_detect(files,'qc')
    files = files[!idx]
    
    for(i in 1:length(files))
    {
        file = files[i]
        obj = readRDS(file)
        a = FetchData(obj,c('nCount_RNA','nFeature_RNA','percent.mt','orig.ident'))
        if(i == 1)
            df = a
        else
            df = rbind(df,a)
    }

    df$orig.ident = factor(df$orig.ident,levels=getDataOrder())

    g = ggplot(df,aes(x=nCount_RNA,y=percent.mt,color=orig.ident)) +
        geom_point(size=.5) +
        facet_wrap(~orig.ident)

    print(g)
    GGSave(plot=g,
           filename='integrationFigures/percent.mt_vs_nCount_RNA.pdf',
           height=12,width=18,units='in')

    h = ggplot(df,aes(x=nFeature_RNA,y=percent.mt,color=orig.ident)) +
        geom_point(size=.5) +
        facet_wrap(~orig.ident)

    dev.new()
    print(h)
    GGSave(plot=h,
           filename='integrationFigures/percent.mt_vs_nFeature_RNA.pdf',
           height=12,width=18,units='in')
}

 
## ####################################################
plotTradeoffs = function()
{
    plotCountTradeoffs()
    plotFeatureTradeoffs()
}
   

## ####################################################
plotCountTradeoffs = function()
{
    files = Sys.glob('unintegrated/*')
    idx = str_detect(files,'qc')
    files = files[!idx]
    
    mtCutoff = 10
    countCutoff = seq(from=0,to=1000,by=100)
    for(i in 1:length(files))
    {
        obj = readRDS(files[i])
        numCells = ncol(obj)

        passRate = c()
        remainingCells = c()
        for(j in 1:length(countCutoff))
        {
            idx = obj$percent.mt <= mtCutoff &
                obj$nCount_RNA >= countCutoff[j]
            passRate[j] = sum(idx) / numCells
            remainingCells[j] = sum(idx)
        }

        a = data.frame(countCutoff,
                       remainingCells,
                       passRate,
                       data=unique(obj$orig.ident))
        if(i == 1)
            df = a
        else
            df = rbind(df,a)
    }
    df$data = factor(df$data,levels=getDataOrder())
                     
    g = ggplot(df,aes(x=countCutoff,y=passRate,color=data,group=data)) +
        geom_line() +
        geom_point() +
        ggtitle('pass rate vs. count cutoff')

    dev.new()
    print(g)
    GGSave(plot=g,
           filename='integrationFigures/passRateVsCountCutoff.pdf',
           height=8,width=10,units='in')

    nameAndMakeDir('integrationTables')
    fileName='integrationTables/percent.mt_vs_nCount_RNA.txt'
    Write.Table(df,fileName)
    
}

## ####################################################
plotFeatureTradeoffs = function()
{
    files = Sys.glob('unintegrated/*')
    idx = str_detect(files,'qc')
    files = files[!idx]
    
    mtCutoff = 10
    featureCutoff = seq(from=0,to=1000,by=100)
    for(i in 1:length(files))
    {
        obj = readRDS(files[i])
        numCells = ncol(obj)

        passRate = c()
        remainingCells = c()
        for(j in 1:length(featureCutoff))
        {
            idx = obj$percent.mt <= mtCutoff &
                obj$nFeature_RNA >= featureCutoff[j]
            passRate[j] = sum(idx) / numCells
            remainingCells[j] = sum(idx)
        }

        a = data.frame(featureCutoff,
                       remainingCells,
                       passRate,
                       data=unique(obj$orig.ident))
        if(i == 1)
            df = a
        else
            df = rbind(df,a)
    }
    df$data = factor(df$data,levels=getDataOrder())
                     
    g = ggplot(df,aes(x=featureCutoff,y=passRate,color=data,group=data)) +
        geom_line() +
        geom_point() +
        ggtitle('pass rate vs. feature cutoff')

    dev.new()
    print(g)
    GGSave(plot=g,
           filename='integrationFigures/passRateVsFeatureCutoff.pdf',
           height=8,width=10,units='in')

    nameAndMakeDir('integrationTables')
    fileName='integrationTables/percent.mt_vs_nFeature_RNA.txt'
    Write.Table(df,fileName)
    
}
             


## ####################################################
getDataOrder = function()
{
    tags = c('2dpf',
             '3dpf',
             'adult19',
             'adult20')
    return(tags)
}



## ####################################################
performQC = function()
{
    files = Sys.glob('unintegrated/*')
    idx = str_detect(files,'qc')
    files = files[!idx]

    mtCutoff = 10

    writeLines(c('mt cutoff 10%',
                 'adult feature cutoff 300',
                 'others 600'),
                 'integrationTables/qc_cutoffs.txt')

    for(file in files)
    {
        obj = readRDS(file)
        tag = unique(obj$orig.ident)
        if(str_detect(tag,'adult'))
            featureCutoff = 300
        else
            featureCutoff = 600

        idx = obj$percent.mt <= mtCutoff &
            obj$nFeature_RNA >= featureCutoff

        obj = obj[,idx]
        fileOut = str_replace(file,'\\.rds','_qc.rds')

        saveRDS(obj,fileOut)
    }
}

## ####################################################
runBasicAnalyses = function(obj)
{
    obj = ScaleData(obj)
    obj = FindVariableFeatures(obj)
    obj = RunPCA(obj)
    obj = RunUMAP(obj,dims=1:20)
    obj = RunTSNE(obj)
    obj = FindNeighbors(obj)

    return(obj)
}

## ####################################################
makeLarvalObject = function()
{
    day2 = readRDS('unintegrated/2dpf_qc.rds')
    day3 = readRDS('unintegrated/3dpf_qc.rds')

    larval = merge(day2,y=day3,
                   project='larval')

    larval = NormalizeData(larval)
    larval = runBasicAnalyses(larval)

    saveDir = nameAndMakeDir('intermediate')
    saveRDS(larval,'intermediate/larval.rds')
}

## ####################################################
makeIntegratedLarvalObject = function()
{
    day2 = readRDS('unintegrated/2dpf_qc.rds')
    day3 = readRDS('unintegrated/3dpf_qc.rds')

    day2 = NormalizeData(day2)
    day2 = FindVariableFeatures(day2)

    day3 = NormalizeData(day3)
    day3 = FindVariableFeatures(day3)

    anchors = FindIntegrationAnchors(list(day2,day3))
    larval = IntegrateData(anchors)

    larval = runBasicAnalyses(larval)  

    saveRDS(larval,
            'intermediate/larval_integrated.rds')
}

## ####################################################
makeAdultObject = function()
{
    adult19 = readRDS('unintegrated/adult19_qc.rds')
    adult20 = readRDS('unintegrated/adult20_qc.rds')

    adult = merge(adult19,y=adult20,
                   project='audlt')

    adult = NormalizeData(adult)
    adult = runBasicAnalyses(adult)

    saveDir = nameAndMakeDir('intermediate')
    saveRDS(adult,'intermediate/adult.rds')
}

## ####################################################
makeCombinedObject = function(integratedLarval,
                              whichKeep,
                              whichAdult='both')
{
    adult = readRDS('intermediate/adult_CellType.rds')

    if(whichAdult == 20)
    {
        idx = adult$orig.ident == 'adult19'
        adult = adult[,idx]
    }

    if(integratedLarval)
    {
        larval = readRDS('intermediate/larval_integrated_CellType.rds')
    } else {
        larval = readRDS('intermediate/larval_CellType.rds')
    }

    
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

    anchors = FindIntegrationAnchors(list(adult,larval))
    combined = IntegrateData(anchors)

    combined = runBasicAnalyses(combined)  

    if(integratedLarval)
        lTag = '_integratedLarval'
    else
        lTag = '_mergedLarval'

    if(whichKeep == 1)
        kTag = '_keepCells1'
    else
        kTag = '_keepCells2'

    fileName = paste0('intermediate/combined',
                      lTag,
                      kTag,
                      '.rds')

    if(whichAdult==20)
            fileName = paste0('intermediate/combined',
                      lTag,
                      kTag,
                      '_adult20.rds')


    saveRDS(combined,fileName)
}

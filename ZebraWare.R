library(Seurat)
library(HandyPack)

## ####################################################
## ####################################################
getSeuratObject = function(which)
{
    fileName = paste0('SeuratObject/',
                      which,
                      '.rds')
    f = readRDS(fileName)

    return(f)
}

## ####################################################
largeGeneList = function()
{
    df = Read.Table('largeGeneList.txt')
    genes = df$gene

    idx = genes == ''
    genes = genes[!idx]

    return(genes)
}

## ####################################################
updateGeneList = function(genesToAdd)
{
    Tic('updating gene list')
    fileName = 'largeGeneList.txt'
    df = Read.Table(fileName)
    Write.Table(df,
                'largeGeneList.save')
    
    gene = c(df$gene,genesToAdd)
    gene = gene[!duplicated(gene)]
    gene = gene[order(gene)]

    df = data.frame(gene)

    Write.Table(df,
                fileName)
    toc()
}

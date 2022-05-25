
library(Seurat)
library(HandyPack)
library(ggplot2)
library(stringr)

rm(list=ls())
graphics.off()

## ####################################################
## ####################################################

## ####################################################
file = '~/TiffanyZebraFish04/unintegrated/adultTogether.rds'
f = readRDS(file)
vars = c('percent.mt','nFeature_RNA','nCount_RNA')

df = FetchData(f,vars)

fStef = readRDS('~/TiffanyZebraFish04/unintegrated/stefan20.rds')
df2 = FetchData(fStef,vars)

df$which = 'together'
df2$which = 'stefan20'

df = rbind(df,df2)

for(var in vars)
{
    title = paste(var,'in rematch and stefan')
    g = ggplot(df,aes_string(x=var)) +
        geom_histogram(bins=50) +
        facet_wrap(~which,ncol=1) +
        ggtitle(title)

    dev.new()
    print(g)

    title = str_replace(title,' ','_')
    saveDir = nameAndMakeDir(c('counts','rawObjects'))
    fileName = paste0(saveDir,'/',title,'.jpg')
    ggsave(plot=g,
           filename=fileName,
           height=8,width=14,units='in')
}

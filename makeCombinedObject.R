
library(Seurat)
library(HandyPack)

rm(list=ls())

source('IntegrationWare.R')
source('CopyWare.R')
## ####################################################
## ####################################################

for(integratedLarval in c(TRUE,FALSE))
    for(whichKeep in 1:2)
    {
        writeLines(paste('integratedLaraval',integratedLarval,whichKeep))
        makeCombinedObject(integratedLarval,whichKeep)
    }






rm(list=ls())
graphics.off()

source('IntegrationWare.R')
source('CopyWare.R')
## ####################################################
## ####################################################
Source = function(x)
{
    source(x)

    ## Make sure these aren't over-written:
    source('IntegrationWare.R')
    source('CopyWare.R')
}
    
makeSeuratObjectFrom10X()

makeQCFigures()

plotTradeoffs()

performQC()

makeLarvalObject()

makeIntegratedLarvalObject()

copyCellTypeIntoNewLarval()

copyCellTypeIntoNewIntegratedLarval()

Source('makeLarvalQCFigures.R')

makeAdultObject()

copyCellTypeIntoNewAdult()

copyCellTypeIntoNewAdult20()

Source('makeAdultQCFigures.R')



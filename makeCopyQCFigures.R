
library(Seurat)
library(HandyPack)
library(plotly)
library(ggplot2)
library(stringr)

rm(list=ls())
graphics.off()

## ####################################################
## ####################################################
files = Sys.glob('intermediateCopyQC/*Cell*')

saveDir = nameAndMakeDir('copyQCFigures')
vars = c('UMAP_1','UMAP_2','orig.ident','CellType')
for(file in files)
{
    f = readRDS(file)
    df = FetchData(f,vars)
    for(var in c('orig.ident','CellType'))
    {
        N = length(unique(df[,var]))
        colors = polychrome(N)
        title = str_replace(file,'intermediateCopyQC/','')
        title = str_replace(title,'.rds','')
        title = paste(var,'in',title)
        g = ggplot(df,aes_string(x='UMAP_1',y='UMAP_2',color=var)) +
            geom_point() +
            scale_color_manual(values=colors) +
            ggtitle(title)

        ## dev.new()
        ## print(g)
        fileName = paste0('copyQCFigures/',
                          title,
                          '.jpg')
        fileName = str_replace_all(fileName,' ','_')
        
        ggsave(plot=g,
               filename=fileName,
               width=14,height=9,units='in')

        p = ggplotly(g)
        print(p)
        fileName = str_replace(fileName,'jpg','html')
        htmlwidgets::saveWidget(as_widget(p), fileName)
    }
}
    

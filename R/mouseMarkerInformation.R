# show markers
# jiang
# 2017/02/17

library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(amap)
library(xlsx)
library(gplots)
library(ggplot2)
library(lattice)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
markerFolder = paste0(adipocyte_dir,'BrW_Markers/')


## biomart get genes
ensemblMartMouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                           host = "dec2016.archive.ensembl.org")


############################
### read necessary dataset
############################
mouseMarker = read.delim(file = paste0(markerFolder,'marker.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Ensembl = rownames(mouseMarker)

## add markers type
mouseMarkerBat = read.delim(file = paste0(markerFolder,'marker.BAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarkerWat = read.delim(file = paste0(markerFolder,'marker.WAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Type = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'BAT','WAT'))

## add markers subcellular location, wether it is mitochondrial or not
# add mitocarta information, downloaded from Mitocarta: https://www.broadinstitute.org/scientific-community/science/programs/metabolic-disease-program/publications/mitocarta/mitocarta-in-0
mitocarta = read.xlsx(file = paste0(markerFolder,'Mouse.MitoCarta2.0.xls'),
                      sheetIndex = 2,stringsAsFactors = FALSE)
mouseMarker$MitoCarta = sapply(toupper(mouseMarker$GeneName),function(x) ifelse(x %in% toupper(mitocarta$Symbol),'Yes','No'))
mouseMarker = mouseMarker[order(mouseMarker$Type,abs(mouseMarker$logFC)),]

write.table(mouseMarker,file = paste0(markerFolder,'mouseMarkerInformation.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)

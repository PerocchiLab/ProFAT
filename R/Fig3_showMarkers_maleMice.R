# confirm markers expression via analyzing male mice data
# jiang
# 2017/07/17

library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(amap)
library(xlsx)
library(gplots)
library(ggplot2)
library(lattice)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
#code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))

## biomart get genes
ensemblMartMouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                           host = "dec2016.archive.ensembl.org")
markerFolder = paste0(adipocyte_dir,'BrW_Markers/')


############################
### read necessary dataset
############################
mouseMarker = read.delim(file = paste0(markerFolder,'marker.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Ensembl = rownames(mouseMarker)

# mouseMarker = merge(mouseMarker,getBM(attributes = c('ensembl_gene_id','external_gene_name'),
#                                       filters = 'ensembl_gene_id',
#                                       values = mouseMarker$Ensembl,
#                                       mart = ensemblMartMouse),
#                     by.x='Ensembl',by.y='ensembl_gene_id')

## add markers type
mouseMarkerBat = read.delim(file = paste0(markerFolder,'marker.BAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarkerWat = read.delim(file = paste0(markerFolder,'marker.WAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Type = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'BAT','WAT'))

## add markers subcellular location
# add mitocarta information
mitocarta = read.xlsx(file = paste0(markerFolder,'Mouse.MitoCarta2.0.xls'),
                      sheetIndex = 2,stringsAsFactors = FALSE)
mouseMarker$MitoCarta = sapply(toupper(mouseMarker$GeneName),function(x) ifelse(x %in% toupper(mitocarta$Symbol),'Yes','No'))
mouseMarker = mouseMarker[order(mouseMarker$Type,abs(mouseMarker$logFC)),]

write.table(mouseMarker,file = paste0(markerFolder,'mouseMarkerInformation.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)


#############################################################################################
### RNAseq data from Martin for marker from male mice
#############################################################################################

# ## female mice
# gender = 'female'
# martinData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m00R_BrBeW/tg/rlog.txt'),
#                         sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
# martinData$Ensembl = rownames(martinData)
# martinDataMarker = merge(mouseMarker,martinData,by='Ensembl')
# colors=c('brown2','brown','papayawhip','orange'),

## male mice
#### boxplot for male data, using different colors, only 2 colors(2 conditions only, BAT and iWAT)
boxplotMarkers = function(markerData,mitoList=mitocarta,cex.lab=0.5,...){
  for(tmp_i in 1:(ncol(markerData)-1)){
    boxplot(markerData[,tmp_i]~markerData[,ncol(markerData)],
            data=markerData,main=colnames(markerData)[tmp_i],
            col=c('brown2','papayawhip'),
            col.main = ifelse(colnames(markerData)[tmp_i] %in% mitoList$Symbol,'blue4','deepskyblue2'),
            ylab='',...)
    title(ylab = "Gene expression (rlog)", line = 2.5, cex.lab = 1.2*cex.lab)
  }
}
# read data
gender = 'male'
maleData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m20R_BrW/tg/rlog.txt'),
                      sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
maleData$Ensembl = rownames(maleData)
martinDataMarker = merge(mouseMarker,maleData,by='Ensembl')
##colors=c('brown2','papayawhip')


# processing merged data
rownames(martinDataMarker) = martinDataMarker$GeneName
martinDataMarker = martinDataMarker[order(martinDataMarker$Type,-(abs(martinDataMarker$logFC))),]
martinDataMarker = as.data.frame(t(martinDataMarker[,12:ncol(martinDataMarker)]))
martinDataMarker$condition  = gsub("\\_R[0-9]+\\(r\\d+\\)$", "", rownames(martinDataMarker))


martinDataMarker4KnownBeBr = martinDataMarker[,colnames(martinDataMarker) %in% 
                                                c('Ucp1','Cidea', 'Cox7a1', 'Pdk4','condition')]
martinDataMarker4UnknownBeBr = martinDataMarker[,colnames(martinDataMarker) %in% 
                                                c('Acaa2', 'Slc25a20', 'Aco2', 'Gm13910','condition')]
martinDataMarker4SpecificBr = martinDataMarker[,colnames(martinDataMarker) %in% 
                                                  c('Zic1', 'Impdh1', 'Tmem246', 'Shmt1','condition')]
martinDataMarker4White = martinDataMarker[,colnames(martinDataMarker) %in% 
                                                 c('Hoxc8', 'Alcam', 'Ar', 'Sgpp1','condition')]
martinDataMarkerSupplementary = martinDataMarker[,!colnames(martinDataMarker) %in% 
                                                   c('Ucp1','Cidea', 'Cox7a1', 'Pdk4',
                                                     'Acaa2', 'Slc25a20', 'Aco2', 'Gm13910',
                                                     'Zic1', 'Impdh1', 'Tmem246', 'Shmt1',
                                                     'Hoxc8', 'Alcam', 'Ar', 'Sgpp1')]
###########################33################
## plot main figure
#############################################
pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationMainFigure_', gender, '.pdf'),
    width = 7,height = 9)
#par(mfrow=c(4,4),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(mfrow=c(4,4), las=1)
boxplotMarkers(martinDataMarker4KnownBeBr,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
boxplotMarkers(martinDataMarker4UnknownBeBr,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
boxplotMarkers(martinDataMarker4SpecificBr,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
boxplotMarkers(martinDataMarker4White,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
dev.off()

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationSupplementary_', gender, '.pdf'),
    width = (183/25.4),height = (247/25.4))
#par(mfrow=c(4,4),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(mfrow=c(4,4), las=1)
boxplotMarkers(martinDataMarkerSupplementary,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
dev.off()

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationSupplementary1page_', gender, '.pdf'),
    width = (183/25.4),height = (247/25.4))
par(mfrow=c(8,6),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
boxplotMarkers(martinDataMarkerSupplementary,
               boxwex=0.5,cex.axis=0.5,cex.main=0.8,lwd=0.2,cex.lab=0.5,
               border = 'gray0',las=2)
dev.off()


pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmation_', gender, '.pdf'),
    width = (183/25.4),height = (247/25.4))
par(mfrow=c(6,5),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
boxplotMarkers(martinDataMarker,
               boxwex=0.5,cex.axis=0.5,cex.main=0.8,lwd=0.2,
               border = 'gray0',las=2)
dev.off()


##############################################################################
#####  plot brown and white markers sepatately
############################################################################
# pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationBAT.pdf'),
#     width = 8,height = 10)
# #par(mfrow=c(4,4),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
# par(mfrow=c(4,4), las=1)
# for(col_i in 1:(ncol(martinDataMarker)-7)){
#   boxplot(martinDataMarker[,col_i]~martinDataMarker[,ncol(martinDataMarker)],
#           data=martinDataMarker,main=colnames(martinDataMarker)[col_i],
#           col=c('brown2','brown','papayawhip','orange'),
#           boxwex=0.5,cex.axis=0.68,cex.main=1,lwd=0.2,border = 'gray0',
#           ylab='Gene expression (rlog)')
# }
# dev.off()
# pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationWAT.pdf'),
#     width = 8,height = 6)
# #par(mfrow=c(2,3),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
# par(mfrow=c(2,3), las=1)
# for(col_i in (ncol(martinDataMarker)-6):(ncol(martinDataMarker)-1)){
#   boxplot(martinDataMarker[,col_i]~martinDataMarker[,ncol(martinDataMarker)],
#           data=martinDataMarker,main=colnames(martinDataMarker)[col_i],
#           col=c('brown2','brown','papayawhip','orange'),
#           boxwex=0.5,cex.axis=0.68,cex.main=1,lwd=0.2,border = 'gray0',
#           ylab='Gene expression (rlog)')
# }
# dev.off()
#############################################################################################

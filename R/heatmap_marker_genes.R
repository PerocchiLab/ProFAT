# show markers
# jiang
# 2017.02.17, 
# modified 2017.11.07

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


### read mouse marker genes list
mouseMarker = read.delim(file = paste0(markerFolder,'mouseMarkerInformation.txt'),
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                         check.names = FALSE)

### get gene expression of markers
mouseExpression = read.delim(file = paste0(markerFolder,'rbatch.data_BrW.txt'),
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                             check.names = FALSE)


## calculate gene expression means within biological replicates
mouseExpressionMean = as.data.frame(t(mouseExpression))
mouseExpressionMean$samples = gsub("\\(r\\d+\\)$", "", rownames(mouseExpressionMean))
mouseExpressionMean = aggregate(.~samples,data = mouseExpressionMean,mean)
rownames(mouseExpressionMean)=mouseExpressionMean$samples
mouseExpressionMean = mouseExpressionMean[,-1] # remove 'sample' column
mouseExpressionMean = as.data.frame(t(mouseExpressionMean))
mouseExpressionMean$Ensembl = rownames(mouseExpressionMean)
mouseMarkerExpressionMean = merge(mouseMarker,mouseExpressionMean,by='Ensembl')
rownames(mouseMarkerExpressionMean) = mouseMarkerExpressionMean$GeneName
mouseMarkerExpressionMean = mouseMarkerExpressionMean[order(mouseMarkerExpressionMean$Type,-abs(mouseMarkerExpressionMean$logFC)),]

## plot the foldchange as barplot
markerFC = mouseMarkerExpressionMean[,c("logFC","Type","GeneName")]
markerFC = markerFC[order(markerFC$Type,-(abs(markerFC$logFC)),decreasing = TRUE),]
markerFC$row = 1:nrow(markerFC)

## plot foldchange
#pdf(file=paste0(adipocyte_dir,'BrW_Markers/fig3markersFoldChange.pdf'),
#    onefile = FALSE,family = "Helvetica",width = 2,height = 8)
barchart(markerFC$row~markerFC$logFC,
         data = markerFC, origin = 0, horizontal = TRUE,
         xlab=list(label = 'log2(FoldChange)',cex=0.7),
         ylab=list(label = 'Genes',cex = 0.7),
         scales=list(y=list(at=markerFC$row,labels=markerFC$GeneName,cex=.6),
                     x=list(cex=0.5))
)
#dev.off()


mouseMarkerExpressionMean = mouseMarkerExpressionMean[,11:ncol(mouseMarkerExpressionMean)]

# pdf(file=paste0(adipocyte_dir,'BrW_Markers/fig3markersHeatmap.pdf'),
#     onefile = FALSE,family = "Helvetica",width = 11,height = 10)
heatmap.2(as.matrix(mouseMarkerExpressionMean[,-1]),
          main = "Markers Expression", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(12,9),     # widens margins around plot
          col = colorRampPalette(c('green3','yellow','red'))(n = 60),  
          #labRow = TRUE,
          keysize = 1.5,srtCol=45,
          cexRow = 1,cexCol = 1,scale = 'row',# if scale = 'row'|'col', the data will be change to Z score instead of the original data
          RowSideColors = unname(sapply(mouseMarkerExpressionMean$MitoCarta,
                                        function(x) ifelse(x == 'Yes','black','gray'))),
          #Rowv = as.dendrogram(row_cluster), # apply default clustering method
          #Colv = as.dendrogram(col_cluster),
          Rowv=NA,dendrogram = 'col')
legend("topright",
       legend=c('Mitocarta','Non-mitocarta'),
       fill=c('black','gray'),
       border=FALSE, bty="n", y.intersp = 1, cex=0.7)
#dev.off()

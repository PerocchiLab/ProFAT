#2017/02/08

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(amap)
library(xlsx)
library(gplots)

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
mouseMarker$mitocarta = sapply(toupper(mouseMarker$GeneName),function(x) ifelse(x %in% toupper(mitocarta$Symbol),'Yes','No'))
mouseMarker = mouseMarker[order(mouseMarker$Type,abs(mouseMarker$logFC)),]

write.table(mouseMarker,file = paste0(markerFolder,'mouseMarkerInformation.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)
########################################################################################
## markers expression 
####################################################################################
mouseExpression = read.delim(file = paste0(markerFolder,'rbatch.data_BrW.txt'),
                             sep = '\t',header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

## get means
mouseExpressionMean = as.data.frame(t(mouseExpression))
mouseExpressionMean$samples = gsub("\\(r\\d+\\)$", "", rownames(mouseExpressionMean))
mouseExpressionMean = aggregate(.~samples,data = mouseExpressionMean,mean)
rownames(mouseExpressionMean)=mouseExpressionMean$samples
mouseExpressionMean = mouseExpressionMean[,-1]
mouseExpressionMean = as.data.frame(t(mouseExpressionMean))
mouseExpressionMean$Ensembl = rownames(mouseExpressionMean)
mouseMarkerExpressionMean = merge(mouseMarker,mouseExpressionMean,by='Ensembl')
rownames(mouseMarkerExpressionMean) = mouseMarkerExpressionMean$GeneName
mouseMarkerExpressionMean = mouseMarkerExpressionMean[order(mouseMarkerExpressionMean$Type,-abs(mouseMarkerExpressionMean$logFC)),]
mouseMarkerExpressionMeanSelected = mouseMarkerExpressionMean[,12:ncol(mouseMarkerExpressionMean)]
mat = apply(mouseMarkerExpressionMeanSelected, 1, function(x) (x-mean(x))/sd(x))
mat=t(mat)
foldChange = mouseMarkerExpressionMean$logFC
location = data.frame(Mitocarta=mouseMarkerExpressionMean[,'mitocarta'])


pdf(file=paste0(adipocyte_dir,'BrW_Markers/fig3aMarkerHeatmap.pdf'),
    onefile = FALSE,family = "Helvetica",width = 11,height = 10)
row_ha = rowAnnotation(foldChange = row_anno_barplot(foldChange, axis = TRUE, 
                                                     axis_side = "bottom",
                                                     axis_gp = gpar(fontsize = 6),
                                                     axis_direction = "reverse",
                                                     axis_label='tt',
                                               gp = gpar(fill = 'lightblue')),
                       width = unit(1.5, "cm"))


ha = rowAnnotation(df = location, 
                   col = list(Mitocarta = c("Yes" = "gray", "No" = "black")),
                   width = unit(0.2, "cm"))

row_ha + ha + Heatmap(mat,cluster_rows=FALSE,row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = colorRampPalette(c('green','yellow','red'))(n =20),  
        heatmap_legend_param = list(title = 'Z-score', color_bar = "continuous"))

dev.off()

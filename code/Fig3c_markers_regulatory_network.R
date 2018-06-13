# regulatory network for markers
# first use iRegulon to get the TFs and targets gene then use R to extract them
# jiang
# 2017/02/24

library(RColorBrewer)
library(biomaRt)
library(xlsx)
library(pheatmap)
library(amap)
library(gplots)
library(ggplot2)
library(lattice)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
#code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))

markerFolder = paste0(adipocyte_dir,'BrW_Markers/')
## biomart get genes

ensemblMartMouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                           host = "dec2016.archive.ensembl.org")
mitocarta = read.xlsx(file = paste0(markerFolder,'Mouse.MitoCarta2.0.xls'),
                      sheetIndex = 2,stringsAsFactors = FALSE)

#####
 getiRegulonRegulatoryNetwork(iregulonFile = paste0(markerFolder,'iregulon/iregulonDefault.tsv'),
                             markerInformationFile = paste0(markerFolder,'mouseMarkerInformation.txt'),
                             logFCcutoff=1,pvalueCutoff = 0.01,
                             sigFile=paste0(markerFolder,'df_exp.txt'),
                             mitoCartaList = mitocarta$Symbol)
 test1=getiRegulonRegulatoryNetwork(iregulonFile = paste0(markerFolder,'iregulon/iregulonDefault.tsv'),
                             markerInformationFile = paste0(markerFolder,'mouseMarkerInformation.txt'),
                             logFCcutoff=1.5,pvalueCutoff = 0.01,
                             sigFile=paste0(markerFolder,'df_exp.txt'),
                             mitoCartaList = mitocarta$Symbol)
getiRegulonRegulatoryNetwork(iregulonFile = paste0(markerFolder,'iregulon/iregulonDefault.tsv'),
                             markerInformationFile = paste0(markerFolder,'mouseMarkerInformation.txt'),
                             logFCcutoff=2,pvalueCutoff = 0.01,
                             sigFile=paste0(markerFolder,'df_exp.txt'),
                             mitoCartaList = mitocarta$Symbol)

test1$selectedTF[toupper(test1$selectedTF) %in% toupper(mouseMarker$GeneName)]


######################################################################################33
##### draw legend via heatmap.2

mouseMarker = read.delim(file = paste0(markerFolder,'mouseMarkerInformation.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
rownames(mouseMarker) = mouseMarker$GeneName

tmpFC = mouseMarker[,c("logFC",'GeneName')]
tmpFCs = cbind(tmpFC[,1],
#               runif(nrow(tmpFC),min = min(tmpFC$logFC),max=max(tmpFC$logFC)),
tmpFC[,1])
## plot foldchange
#pdf(file=paste0(adipocyte_dir,'iregulon/legend1.pdf'),
pdf(file=paste0(plots_dir,'legend3.pdf'),
    onefile = FALSE,family = "Helvetica",width = 8,height = 8)
heatmap.2(t(as.matrix(tmpFCs)),
          main = "Network Legend", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          #margins = c(12,9),     # widens margins around plot
          col = colorRampPalette(c('green2','white','red'))(n = 20),       # use on color palette defined earlier
          keysize = 4,srtCol=45,symbreaks = TRUE,symm = FALSE,
          cexRow = 1,cexCol = 1,
          Rowv=NA,dendrogram = 'col')
# pheatmap(as.matrix(rbind(tmpFC$logFC,tmpFC$logFC)),
#          col = colorRampPalette(c('green2','white','red'))(n = 20))

dev.off()

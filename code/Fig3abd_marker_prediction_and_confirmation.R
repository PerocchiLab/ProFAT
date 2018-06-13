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

# #####################################################################
# ### add uniprot subcellular location
# #####################################################################
# ## first using python to download localization
# for(gene_i in 1:nrow(mouseMarker)){
#   tmpcommand = paste0('python ',paste0(adipocyte_dir,'MarkerPrediction/query_uniprot.py '),
#                       mouseMarker$Ensembl[gene_i],' ', paste0(adipocyte_dir,'MarkerPrediction/location/'),
#                       paste0(' ',mouseMarker$Ensembl[gene_i],'.txt'))
#   system(tmpcommand)
# }
# 
# mouseMarker$Location = sapply(mouseMarker$Ensembl,function(x) get_subcellular_location(x))
# sum(mouseMarker$Location !="") ##55, 4 do not have subcellular location information
# mouseMarker$Mitochondria = as.integer(grepl('Mitochondrion',mouseMarker$Location))
# mouseMarker$Cytoplasm = as.integer(grepl('Cytoplasm',mouseMarker$Location))
# mouseMarker$Nucleus = as.integer(grepl('Nucleus',mouseMarker$Location))
# mouseMarker$ER = as.integer(grepl('Endoplasmic',mouseMarker$Location))
# write.table(mouseMarker,file = paste0(markerFolder,'mouseMarkerInformation.txt'),
#             sep = '\t',row.names = FALSE,quote = FALSE)

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

markerFC = mouseMarkerExpressionMean[,c("logFC","Type","GeneName")]
markerFC = markerFC[order(markerFC$Type,-(abs(markerFC$logFC)),decreasing = TRUE),]
markerFC$row = 1:nrow(markerFC)
## plot foldchange
pdf(file=paste0(adipocyte_dir,'BrW_Markers/fig3markersFoldChange.pdf'),
    onefile = FALSE,family = "Helvetica",width = 2,height = 8)
barchart(markerFC$row~markerFC$logFC,
         data = markerFC, origin = 0, horizontal = TRUE,
         xlab=list(label = 'log2(FoldChange)',cex=0.7),
         ylab=list(label = 'Genes',cex = 0.7),
         scales=list(y=list(at=markerFC$row,labels=markerFC$GeneName,cex=.6),
                     x=list(cex=0.5))
         )
dev.off()


mouseMarkerExpressionMean = mouseMarkerExpressionMean[,11:ncol(mouseMarkerExpressionMean)]

pdf(file=paste0(adipocyte_dir,'BrW_Markers/fig3markersHeatmap.pdf'),
    onefile = FALSE,family = "Helvetica",width = 11,height = 10)
heatmap.2(as.matrix(mouseMarkerExpressionMean[,-1]),
          main = "Markers Expression", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(12,9),     # widens margins around plot
          # col = colorRampPalette(c('lightyellow2','blue'))(n = 50),       # use on color palette defined earlier
          #col = colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160),       # use on color palette defined earlier
          #col = colorRampPalette(c('firebrick2','white','royalblue'))(n = 160),
          col = colorRampPalette(c('green3','yellow','red'))(n = 60),  
         #col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
          #col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
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
dev.off()

###
# plot gene expression without average
mouseExpression$Ensembl = rownames(mouseExpression)
mouseMarkerExpression = merge(mouseMarker,mouseExpression,by='Ensembl')
rownames(mouseMarkerExpression) = mouseMarkerExpression$GeneName
mouseMarkerExpression = mouseMarkerExpression[order(mouseMarkerExpression$Type,-(abs(mouseMarkerExpression$logFC))),]
mouseMarkerExpression = mouseMarkerExpression[,11:ncol(mouseMarkerExpression)]

pdf(file=paste0(adipocyte_dir,'BrW_Markers/markersExpression.pdf'),
    onefile = FALSE,family = "Helvetica",width = 10,height = 10)
heatmap.2(as.matrix(mouseMarkerExpression[,-1]),
          main = "Markers Expression", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins = c(12,9),     # widens margins around plot
          # col = colorRampPalette(c('lightyellow2','blue'))(n = 50),       # use on color palette defined earlier
          col = colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160),       # use on color palette defined earlier
          #col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
          #col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
          #labRow = TRUE,
          keysize = 1,
          cexRow = 1,cexCol = 0.5,scale = 'row',
          RowSideColors = unname(sapply(mouseMarkerExpression$mitocarta,
                                        function(x) ifelse(x == 'Yes','gray','black'))),
          #Rowv = as.dendrogram(row_cluster), # apply default clustering method
          #Colv = as.dendrogram(col_cluster),
          Rowv=NULL,dendrogram = 'col')
legend("topright",
       legend=c('Mitochondria','non-mitochondria'),
       fill=c('gray','black'),
       border=FALSE, bty="n", y.intersp = 1, cex=0.7)
dev.off()

#############################################################################################
### RNAseq data from Martin for marker
#############################################################################################
martinData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m00R_BrBeW/tg/rlog.txt'),
                        sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
martinData$Ensembl = rownames(martinData)
martinDataMarker = merge(mouseMarker,martinData,by='Ensembl')
rownames(martinDataMarker) = martinDataMarker$GeneName
martinDataMarker = martinDataMarker[order(martinDataMarker$Type,-(abs(martinDataMarker$logFC))),]
martinDataMarker = as.data.frame(t(martinDataMarker[,12:ncol(martinDataMarker)]))
martinDataMarker$condition  = gsub("\\_R[0-9]\\(r\\d+\\)$", "", rownames(martinDataMarker))
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
pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationMainFigure.pdf'),
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

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationSupplementary.pdf'),
    width = (183/25.4),height = (247/25.4))
#par(mfrow=c(4,4),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(mfrow=c(4,4), las=1)
boxplotMarkers(martinDataMarkerSupplementary,
               boxwex=0.5,cex.axis=0.55,cex.main=1,lwd=0.2,border = 'gray0')
dev.off()

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markersConfirmationSupplementary1page.pdf'),
    width = (183/25.4),height = (247/25.4))
par(mfrow=c(8,6),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
boxplotMarkers(martinDataMarkerSupplementary,
               boxwex=0.5,cex.axis=0.5,cex.main=0.8,lwd=0.2,cex.lab=0.5,
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
###########################################################################################
## David NCBI
#######################
## GO analysis 
#######################
## get non-redundant GO terms
redudantFiles = list.files(path = paste0(markerFolder,'david/'),recursive = FALSE,
                           full.names = TRUE,pattern = '.txt')
redudantFiles = redudantFiles[grepl('Classification',redudantFiles)]
redudantFiles = redudantFiles[!grepl('Unique',redudantFiles)]
for(file_i in 1:length(redudantFiles)){
  remove_redundant_terms(inFilePath = redudantFiles[file_i])
  print(file_i)
}
 
## plot GO for all clusters seperately
clusterFiles = list.files(path = paste0(markerFolder,'david/'),recursive = FALSE,
                          full.names = TRUE,pattern = '.txt')
clusterFiles = clusterFiles[grepl('Classification',clusterFiles)]
clusterFiles = clusterFiles[grepl('Unique',clusterFiles)]
for(cluster_i in 1:length(clusterFiles)){
  if(grepl('Kegg',clusterFiles[cluster_i])){
    plot_david_output(inFilePath = clusterFiles[cluster_i],
                      figSavePath = gsub('txt$','pdf',clusterFiles[cluster_i]),
                      splitSymbol = ':',top=5,label.cex = 0.6,
                      col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
                      main=basename(clusterFiles[cluster_i]),cex.names = 0.65)
  }else{
    plot_david_output(inFilePath = clusterFiles[cluster_i],
                      figSavePath = gsub('txt$','pdf',clusterFiles[cluster_i]),
                      splitSymbol = '~',top=5,label.cex = 0.6,
                      col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
                      main=basename(clusterFiles[cluster_i]),cex.names = 0.65)
    
  }
}

### plot GO terms on one single pdf ## using the High classification standard
clusterFilesHigh = clusterFiles[grepl('High',clusterFiles)]
clusterFilesHighGO = clusterFilesHigh[!grepl('Kegg',clusterFilesHigh)]
pdf(file=paste0(adipocyte_dir,'BrW_Markers/david/batDavidGO.pdf'),
    onefile = FALSE,family = "Helvetica",width = 12,height = 4)
par(mfrow=c(1,3))
for(cluster_i in 1:length(clusterFilesHighGO)){
  plot_david_output(inFilePath = clusterFilesHighGO[cluster_i],
                    figSavePath = NULL,
                    splitSymbol = '~',top=5,label.cex = 0.65,
                    col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
                    main=ifelse(grepl('BP',clusterFilesHighGO[cluster_i]),'Biological Process',
                                ifelse(grepl('MF',clusterFilesHighGO[cluster_i]),
                                       'Molecular Function','Cellular Component')),
                    cex.names = 0.5)
}
dev.off()

#########################
## pathway kegg david
########################
clusterFilesMedium = clusterFiles[grepl('Medium',clusterFiles)]
clusterFilesMediumKegg = clusterFilesMedium[grepl('Kegg',clusterFilesMedium)]
pdf(file=paste0(adipocyte_dir,'BrW_Markers/david/batDavidPathway.pdf'),
    onefile = FALSE,family = "Helvetica",width = 3,height = 3)
for(cluster_i in 1:length(clusterFilesMediumKegg)){
  plot_david_output(inFilePath = clusterFilesMediumKegg[cluster_i],
                    figSavePath = NULL,
                    splitSymbol = ':',top=5,label.cex = 0.6,
                    col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
                    main='KEGG',
                    cex.names = 0.5)
}
dev.off()

###################################################################################
## cpdb http://cpdb.molgen.mpg.de/MCPDB/ 
##########################################
##  GO from cpdb
##################

batGOcpdb = read.delim(file = paste0(markerFolder,'cpdb/batGo.tab'),
                       sep = '\t',header = TRUE,stringsAsFactors = FALSE)
batGOcpdb = batGOcpdb[batGOcpdb$p.value<0.01,]
termLevels = unique(batGOcpdb$term_level)
termCategory = unique(batGOcpdb$term_category)
# for(level_i in 1:length(termLevels)){
# draw 3 pictures in 3 files
#   batGOcpdbLevelSelected = batGOcpdb[batGOcpdb$term_level==termLevels[level_i],]
#   pdf(file=paste0(markerFolder,'cpdb/batCpdbGOLevel',termLevels[level_i],'.pdf'),
#       onefile = FALSE,family = "Helvetica",width = 10,height =20)
# draw 3 picture in 1 file
  pdf(file=paste0(markerFolder,'cpdb/batCpdbGOLevels.pdf'),
      onefile = FALSE,family = "Helvetica",width = 18,height =12)
  par(mfrow=c(3,3))
  for(level_i in 1:length(termLevels)){
  batGOcpdbLevelSelected = batGOcpdb[batGOcpdb$term_level==termLevels[level_i],]
  for(category_i in 1:length(termCategory)){
    batGOcpdbLevelSelectedCat = batGOcpdbLevelSelected[batGOcpdbLevelSelected$term_category == termCategory[category_i],]
    if(nrow(batGOcpdbLevelSelectedCat)>5){
      batGOcpdbLevelSelectedCat = batGOcpdbLevelSelectedCat[order(batGOcpdbLevelSelectedCat$p.value),]
      batGOcpdbLevelSelectedCat = batGOcpdbLevelSelectedCat[1:5,]
    }
    batGOcpdbLevelSelectedCat = batGOcpdbLevelSelectedCat[order(-batGOcpdbLevelSelectedCat$p.value),]
    tmpbar = barplot(-log(batGOcpdbLevelSelectedCat$p.value),
            horiz = TRUE,
            #names.arg = wrap.label(batGOcpdbLevelSelectedCat$term_name,10),
            beside = TRUE, las=1,
            border = NA,
            xlab = "-log(p-value)",cex.names = 0.55,
            col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
            main=paste0(ifelse(termCategory[category_i]=='b','Biological Process',
                        ifelse(termCategory[category_i]=='m',
                        'Molecular Function','Cellular Component')),' (Level ',termLevels[level_i],')'))
    text(x = -0.05, y = tmpbar, srt = 0,adj = 1, 
         labels = wrap.label(batGOcpdbLevelSelectedCat$term_name,15), 
         xpd = TRUE,cex = 0.65)
  }

}
  dev.off()
  rm(tmpbar)
  #####################################################################
  ## wat cpdb
  watGOcpdb = read.delim(file = paste0(markerFolder,'cpdb/watGo.tab'),
                         sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  watGOcpdb = watGOcpdb[watGOcpdb$p.value<0.01,]
  termLevels = unique(watGOcpdb$term_level)
  termCategory = unique(watGOcpdb$term_category)
  # for(level_i in 1:length(termLevels)){
  # draw 3 pictures in 3 files
  #   batGOcpdbLevelSelected = batGOcpdb[batGOcpdb$term_level==termLevels[level_i],]
  #   pdf(file=paste0(markerFolder,'cpdb/batCpdbGOLevel',termLevels[level_i],'.pdf'),
  #       onefile = FALSE,family = "Helvetica",width = 10,height =20)
  # draw 3 picture in 1 file
  pdf(file=paste0(markerFolder,'cpdb/watCpdbGOLevels.pdf'),
      onefile = FALSE,family = "Helvetica",width = 10,height =4)
  par(mfrow=c(1,2))
  for(level_i in 1:length(termLevels)){
    watGOcpdbLevelSelected = watGOcpdb[watGOcpdb$term_level==termLevels[level_i],]
    for(category_i in 1:length(termCategory)){
      watGOcpdbLevelSelectedCat = watGOcpdbLevelSelected[watGOcpdbLevelSelected$term_category == termCategory[category_i],]
      if(nrow(watGOcpdbLevelSelectedCat)>5){
        watGOcpdbLevelSelectedCat = watGOcpdbLevelSelectedCat[order(watGOcpdbLevelSelectedCat$p.value),]
        watGOcpdbLevelSelectedCat = watGOcpdbLevelSelectedCat[1:5,]
      }
      watGOcpdbLevelSelectedCat = watGOcpdbLevelSelectedCat[order(-watGOcpdbLevelSelectedCat$p.value),]
      tmpbar = barplot(-log(watGOcpdbLevelSelectedCat$p.value),
              horiz = TRUE,
              #names.arg = wrap.label(watGOcpdbLevelSelectedCat$term_name,5),
              beside = TRUE, las=1,
              border = NA,
              xlab = "-log(p-value)",cex.names = 0.55,
              col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[40],
              main=paste0(ifelse(termCategory[category_i]=='b','Biological Process',
                                 ifelse(termCategory[category_i]=='m',
                                        'Molecular Function','Cellular Component')),
                          ' (Level ',termLevels[level_i],')'))
      text(x = -0.05, y = tmpbar, srt = 0,adj = 1, 
           labels = wrap.label(watGOcpdbLevelSelectedCat$term_name,15), 
           xpd = TRUE,cex = 0.65)
    }
    
  }
  dev.off()
  
######################
# cpdb pathway enrichment
#####################
  batPathwayCpdb = read.delim(file = paste0(markerFolder,'cpdb/batPathways.tab'),
                         sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  batPathwayCpdb = batPathwayCpdb[batPathwayCpdb$p.value<0.01,]
  pathwaySource = unique(batPathwayCpdb$source)

  pdf(file=paste0(markerFolder,'cpdb/batCpdbPathways.pdf'),
      onefile = FALSE,family = "Helvetica",width = 12,height =14)
  par(mfrow=c(2,2))
  for(source_i in 1:length(pathwaySource)){
#       pdf(file=paste0(markerFolder,'cpdb/batCpdbPathways',pathwaySource[source_i],'.pdf'),
#           onefile = FALSE,family = "Helvetica",width = 6,height =8)
    batPathwayCpdbSource = batPathwayCpdb[batPathwayCpdb$source == pathwaySource[source_i],]
    if(nrow(batPathwayCpdbSource)>5){
      batPathwayCpdbSource = batPathwayCpdbSource[order(batPathwayCpdbSource$p.value),]
      batPathwayCpdbSource = batPathwayCpdbSource[1:5,]
    }
    batPathwayCpdbSource = batPathwayCpdbSource[order(-batPathwayCpdbSource$p.value),]
   tmpbar= barplot(-log(batPathwayCpdbSource$p.value),
            horiz = TRUE,
            #names.arg = wrap.label(batPathwayCpdbSource$pathway,5),
            beside = TRUE, las=1,
            border = NA,
            xlab = "-log(p-value)",cex.names = 0.5,
            col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
            main=paste0(pathwaySource[source_i],' Pathway'))
   text(x = -0.05, y = tmpbar, srt = 0,adj = 1, 
        labels = wrap.label(batPathwayCpdbSource$pathway,15), 
        xpd = TRUE,cex = 0.65)
    # dev.off()
    }
dev.off()
  rm(tmpbar)
  
######################
  # cpdb pathway enrichment
#####################
batComplexCpdb = read.delim(file = paste0(markerFolder,'cpdb/batComplex.tab'),
                            sep = '\t',header = TRUE,stringsAsFactors = FALSE)
batComplexCpdb = batComplexCpdb[batComplexCpdb$p.value<0.01,]
complexSource = unique(batComplexCpdb$source)
  
    pdf(file=paste0(markerFolder,'cpdb/batCpdbComplex.pdf'),
        onefile = FALSE,family = "Helvetica",width = 14,height =7)
    par(mfrow=c(1,2))
  for(source_i in 1:length(complexSource)){
#     pdf(file=paste0(markerFolder,'cpdb/batCpdbComplex',complexSource[source_i],'.pdf'),
#         onefile = FALSE,family = "Helvetica",width = 9,height =8)
    batComplexCpdbSource = batComplexCpdb[batComplexCpdb$source == complexSource[source_i],]
    if(nrow(batComplexCpdbSource)>5){
      batComplexCpdbSource = batComplexCpdbSource[order(batComplexCpdbSource$p.value),]
      batComplexCpdbSource = batComplexCpdbSource[1:5,]
    }
    batComplexCpdbSource = batComplexCpdbSource[order(-batComplexCpdbSource$p.value),]
   tmpbar= barplot(-log(batComplexCpdbSource$p.value),
            horiz = TRUE,
            #names.arg = wrap.label(batComplexCpdbSource$complex_name,10),
            beside = TRUE, las=1,
            border = NA,
            xlab = "-log(p-value)",cex.names = 0.5,
            col=colorRampPalette(c('firebrick2','lightyellow','royalblue'))(n = 160)[120],
            main=paste0(complexSource[source_i],' Complex'))
     text(x = -0.05, y = tmpbar, srt = 0,adj = 1, 
          labels = wrap.label(batComplexCpdbSource$complex_name,10), 
          xpd = TRUE,cex = 0.5)
    #dev.off()
  }
  dev.off()
rm(tmpbar)

#############################################################################################
### marker RNA-seq evidence
#############################################################################################
martinData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m00R_BrBeW/tg/rlog.txt'),
                        sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
martinData$Ensembl = rownames(martinData)


mouseExpression$Ensembl = rownames(mouseExpression)
mouseMarkerExpression = merge(mouseMarker,mouseExpression,by='Ensembl')
rownames(mouseMarkerExpression) = mouseMarkerExpression$GeneName
mouseMarkerExpression = mouseMarkerExpression[order(mouseMarkerExpression$Type,-(abs(mouseMarkerExpression$logFC))),]
mouseMarkerExpression = mouseMarkerExpression[,11:ncol(mouseMarkerExpression)]
#############################################################################################
### Subcellular Location
#############################################################################################
### using python to download localization
for(gene_i in 1:nrow(mouseMarker)){
  tmpcommand = paste0('python ',paste0(adipocyte_dir,'MarkerPrediction/query_uniprot.py '),
                      mouseMarker$Ensembl[gene_i],' ', paste0(adipocyte_dir,'MarkerPrediction/location/'),
                      paste0(' ',mouseMarker$Ensembl[gene_i],'.txt'))
  system(tmpcommand)
}

#########################
### add subcelluar location to mouseMarker

mouseMarker$Location = sapply(mouseMarker$Ensembl,function(x) get_subcellular_location(x))
sum(mouseMarker$Location !="") ##55, 4 do not have subcellular location information
mouseMarker$Mitochondria = as.integer(grepl('Mitochondrion',mouseMarker$Location))
mouseMarker$Cytoplasm = as.integer(grepl('Cytoplasm',mouseMarker$Location))
mouseMarker$Nucleus = as.integer(grepl('Nucleus',mouseMarker$Location))
mouseMarker$ER = as.integer(grepl('Endoplasmic',mouseMarker$Location))


write.table(mouseMarker,file = paste0(adipocyte_dir,'BrW_Markers/markersSubcellularLocation.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)


write.xlsx(mouseMarker$external_gene_name[grepl('Mitochondrion',mouseMarker$Location)],
           file = paste0(adipocyte_dir,'BrW_Markers/subcellularLocation.xls'),
           append = FALSE,sheetName = 'Mitochondrion',row.names = FALSE,col.names = FALSE)


write.xlsx(mouseMarker$external_gene_name[grepl('Cytoplasm',mouseMarker$Location)],
           file = paste0(adipocyte_dir,'BrW_Markers/subcellularLocation.xls'),
           append = TRUE,sheetName = 'Cytoplasm',row.names = FALSE,col.names = FALSE)

write.xlsx(mouseMarker$external_gene_name[grepl('Nucleus',mouseMarker$Location)],
           file = paste0(adipocyte_dir,'BrW_Markers/subcellularLocation.xls'),
           append = TRUE,sheetName = 'Nucleus',row.names = FALSE,col.names = FALSE)

write.xlsx(mouseMarker$external_gene_name[grepl('Endoplasmic',mouseMarker$Location)],
           file = paste0(adipocyte_dir,'BrW_Markers/subcellularLocation.xls'),
           append = TRUE,sheetName = 'ER',row.names = FALSE,col.names = FALSE)
### new list do not contain proteins locating at Golgi and Lysosome
# 
# write.xlsx(mouseMarker$external_gene_name[grepl('Golgi',mouseMarker$Location)],
#            file = paste0(adipocyte_dir,'MarkerPrediction/subcellularLocation.xls'),
#            append = TRUE,sheetName = 'Golgi',row.names = FALSE,col.names = FALSE)
# 
# write.xlsx(mouseMarker$external_gene_name[grepl('Lysosome',mouseMarker$Location)],
#            file = paste0(adipocyte_dir,'MarkerPrediction/subcellularLocation.xls'),
#            append = TRUE,sheetName = 'Lysosome',row.names = FALSE,col.names = FALSE)

#############################################################################################

##########################################################################
### PPI

markerIDs = getBM(attributes = c('ensembl_gene_id','external_gene_name','ensembl_peptide_id'),
                  filters = 'ensembl_gene_id',
                  values = mouseMarker$Ensembl,
                  mart = ensemblMartMouse)
write.table(markerIDs,file = paste0(adipocyte_dir,'BrW_Markers/markerIDs.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)

mousePPI = read.delim(paste0(adipocyte_dir,'MarkerPrediction/ppi/10090.protein.links.v10.txt'),
                      sep = ' ',header = TRUE,stringsAsFactors = FALSE)  
mousePPI = mousePPI[mousePPI$combined_score>=700,]

mousePPI$markerNumber = colSums(apply(mousePPI,1,function(x) x %in% paste0('10090.',unique(markerIDs$ensembl_peptide_id))))
markerPPI = mousePPI[mousePPI$markerNumber==2,]
markerPPI$ensemblGene1 = ""
markerPPI$ensemblGene2 = ""
markerPPI$geneName1 = ""
markerPPI$geneName2 = ""
for(ppi_i in 1:nrow(markerPPI)){
  markerPPI$ensemblGene1[ppi_i] = from_ensemblPeptide_to_genes(markersMapping = markerIDs,
                                                               ensemblPeptideID = markerPPI$protein1[ppi_i],
                                                               colName = 'ensembl_gene_id')
  markerPPI$ensemblGene2[ppi_i] = from_ensemblPeptide_to_genes(markersMapping = markerIDs,
                                                               ensemblPeptideID = markerPPI$protein2[ppi_i],
                                                               colName = 'ensembl_gene_id')
  markerPPI$geneName1[ppi_i] = from_ensemblPeptide_to_genes(markersMapping = markerIDs,
                                                            ensemblPeptideID = markerPPI$protein1[ppi_i],
                                                            colName = 'external_gene_name')
  markerPPI$geneName2[ppi_i] = from_ensemblPeptide_to_genes(markersMapping = markerIDs,
                                                            ensemblPeptideID = markerPPI$protein2[ppi_i],
                                                            colName = 'external_gene_name')
}

write.table(markerPPI,file = paste0(adipocyte_dir,'BrW_Markers/markersPPI.txt'),
            sep = '\t',row.names = FALSE,quote = FALSE)

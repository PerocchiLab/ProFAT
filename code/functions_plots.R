create_new_dir <- function(dir_path,dir_name=NULL){
  # give the dir path to a variable, if the path does not exist, create it
  new_dir= ifelse(is.null(dir_name),dir_path,paste0(dir_path,dir_name))

  if (!file.exists(new_dir)){
    dir.create(new_dir,recursive=TRUE)
  }
  return (new_dir)
}


read_files_by_line <- function(file_path){
    temp_conn = file(file_path,open="r")
    all_lines = readLines(temp_conn)
    close(temp_conn)
    return(all_lines)
}


#' for a given matrix with gene expression for all samples in 1 studies, caculate the
#' distance among samples first and then display it in a heatmap
#' @param geneExpressionFile a file path contain the gene expression, with each row
#' a gene and each column a sample
#' @param saveFile a file path where to save the heatmap, pdf format
plots_heatmap_distance <- function(geneExpressionFile,sampleIDmap = NULL, saveFile=NULL, dist_metric='euclidean', ...){
  tmpData =  read.delim(file=geneExpressionFile,
                        sep = '\t',header = TRUE,row.names = 1,check.names = FALSE)
  
  colors = colorRampPalette( rev(brewer.pal(8, "Greens")) )(255)
  #colors = colorRampPalette(c('green','orange','red'))(n = 20)
  #colors = colorRampPalette(c("lightgreen", "orange", "red"))(50)
  #colors = colorRampPalette( rev(brewer.pal(8, "Greys")) )(255)
  font_size = ceiling(300/dim(tmpData)[2])
  if (font_size>10){
    font_size = 10
  }
  if(dist_metric %in% c("showGeneValue","showGeneValueScaleColumn")){
    tmpData = as.matrix(t(tmpData))
    # rowname original format 'BAT.line2_hR2'
    rownames(tmpData)= sapply(rownames(tmpData),
                              function(x) ifelse(grepl('_hM|_hR',x),gsub('\\.','\\_',x),x))
    colnames(tmpData)=NULL
    
    if(grepl('ScaleColumn',dist_metric)){
      if(is.null(saveFile)){
        pheatmap(tmpData,
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize_col = font_size,scale = 'column',
                 col=colors, fontsize_row=font_size, ...)
      }else{
        pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
        pheatmap(tmpData,
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize_col = font_size,scale = 'column',
                 col=colors, fontsize_row=font_size, ...)
        dev.off()
      }
    }else{
      if(is.null(saveFile)){
        pheatmap(tmpData,
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize_col = font_size,scale = 'none',
                 col=colors, fontsize_row=font_size, ...)
      }else{
        pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
        pheatmap(tmpData,
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize_col = font_size,scale = 'none',
                 col=colors, fontsize_row=font_size, ...)
        dev.off()
      }
    }
    
    
  }else{
    #datDist  = dist(t(tmpData),method = "euclidean")
    datDist = Dist(t(tmpData), method=dist_metric)
    datDistMatrix = as.matrix(datDist)
    # rowname original format 'BAT.line2_hR2'
    rownames(datDistMatrix)= sapply(rownames(datDistMatrix),
                                    function(x) ifelse(grepl('_hM|_hR',x),gsub('\\.','\\_',x),x))
    colnames(datDistMatrix)=NULL
    if(is.null(saveFile)){
      pheatmap(datDistMatrix,
               clustering_distance_rows=datDist,
               clustering_distance_cols=datDist,
               fontsize_col = font_size,
               col=colors, fontsize_row=font_size, ...)
    }else{
      pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
      pheatmap(datDistMatrix,
               clustering_distance_rows=datDist,
               clustering_distance_cols=datDist,
               fontsize_col = font_size,
               col=colors, fontsize_row=font_size, ...)
      dev.off()
    }
  }
}

get_adipocyte_type = function(sampleName){
  if(grepl('BAT\\(T\\)',sampleName)){
    return('BAT(T)')
  }else if(grepl('BAT',sampleName)){
    return('BAT')
  }else if(grepl('WAT\\(T',sampleName)){
    return('WAT(T)')
  }else{
    return('WAT')
  }
}



get_adipocyte_colors = function(adipocyte_character){
  adipocyte_character = gsub('BAT\\(T\\)','brown',adipocyte_character)
  adipocyte_character = gsub('BAT','brown2',adipocyte_character)
  adipocyte_character = gsub('WAT\\(T\\)','orange',adipocyte_character)
  adipocyte_character = gsub('WAT','papayawhip',adipocyte_character)
  return(adipocyte_character)
}
#' for a given matrix with gene expression, rows are genes and cols are samples
#' we caculate the PCA and display
#' @param geneExpressionFile a file path contain the gene expression, with each row
#' a gene and each column a sample
#' @param saveFile a file path where to save the heatmap, pdf format
#' @param colorGroups is a fator that contain the color for each sample,
#' the order should be consistant with the cols in geneExpressionFile
#' @param pchGroups is a fator that contain the pch  types for each sample,
#' the order should be consistant with the cols in geneExpressionFile
plots_pca <- function(geneExpressionFile,saveFile,colorGroups=NULL, 
                      pchGroups=NULL,ypc=2, ...){
  tmpData =  read.delim(file=geneExpressionFile,sep = '\t',header = TRUE,
                        row.names = 1,check.names = FALSE)
  # if the sample group information is provided by the author,we use the given informaton,
  # otherwise we extract the information from the colnames of the given dataset
  if(is.null(colorGroups)){
    tmpColorGroups = unname(sapply(colnames(tmpData),function(x) get_adipocyte_type(x)))
  }else{
    tmpColorGroups = colorGroups
  }
  if(is.null(pchGroups)){
    #tmpPchGroups = unname(sapply(colnames(tmpData),function(x) strsplit(x,'\\_|\\(r')[[1]][2]))
    tmpPchGroups = str_extract(colnames(tmpData), pattern = '[M|R][0-9]+')
  }else{
    tmpPchGroups = pchGroups
  }
  if(length(unique(tmpPchGroups))==1){
    preDefinePch = c(1)
  }else if(length(unique(tmpPchGroups))>7){
    preDefinePch = rev(c(0:10,12:13))
  }else{
    preDefinePch = c(0:10,12:13)
  }
  tmpDataPCA=prcomp(t(tmpData))
  tmpDataPCAsummary = summary(tmpDataPCA)
  ###add sample type to tmpDataPCA$x
  dfCircle = as.data.frame(tmpDataPCA$x)
  dfCircle$type = as.factor(sapply(row.names(tmpDataPCA$x),function(x) get_adipocyte_type(x)))
  pdf(file = saveFile,onefile = FALSE)
  plot(0, 0, type="n", ann=FALSE, axes=FALSE)
  u <- par("usr") # The coordinates of the plot area
  rect(u[1], u[3], u[2], u[4], border=NA,
       #col=colorRampPalette( rev(brewer.pal(8, "Greens")) )(255)[250]
       col='grey')
  par(new=TRUE)
  plot(tmpDataPCA$x[,1],tmpDataPCA$x[,ypc],
       xlim = c(min(tmpDataPCA$x[,1]-2),max(tmpDataPCA$x[,1]+15)),
       col=get_adipocyte_colors(tmpColorGroups),
       pch=preDefinePch[as.fumeric(tmpPchGroups)],
       xlab = paste0('PC1: ',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC1'],2)*100,'% variance'),
       ylab = paste0('PC',ypc,': ',round(tmpDataPCAsummary$importance["Proportion of Variance",][paste0('PC',ypc)],2)*100,'% variance'),
       ...)
  #   # add legend by groups, on the topleft of the figure
#   legend('bottom',unique(tmpColorGroups),col = get_adipocyte_colors(unique(tmpColorGroups)),
  legend('bottom',unique(tmpColorGroups),col = get_adipocyte_colors(unique(tmpColorGroups)),
         lty = 1,bty = 'n',text.col = get_adipocyte_colors(unique(tmpColorGroups)))
  legend('right',unique(tmpPchGroups),col = 'black',bty = 'n',
         pch = unique(preDefinePch[as.fumeric(tmpPchGroups)]),
         cex = 1)
  #     for(tmp_i in 1:nrow(tmpDataPCA$x)){
  #       text(tmpDataPCA$x[tmp_i,1],(tmpDataPCA$x[tmp_i,ypc]+2),rownames(tmpDataPCA$x)[tmp_i],
  #            col=get_adipocyte_colors(tmpColorGroups)[tmp_i],cex=0.3)
  #     }
  
  #
  
#   ## add circles
#   with(dfCircle, dataEllipse(PC1, PC2, type, 
#                              #labels=rownames(dfCircle),
#                              draw = TRUE,add = TRUE,plot.points=FALSE,
#                              col = c('brown','blue','green'),lwd=0.5,
#                              fill=TRUE,fill.alpha = 0.1, # if transparency and degree
#                              center.pch=15,radius=radius,
#                              levels = 0.99,robust = TRUE, # use cov.trob function calculate the center and covariance matrix
#                              group.labels=NULL))
  
  dev.off()
}


#################################################3
## functions in plot_marker.R
##################################################
check_if_mouse_markers = function(testGeneNames,markerList){
    if(testGeneNames==''|is.na(testGeneNames)){
        return(FALSE)
    }else{
        #testGeneNames = 'ICL At3g21720 MSD21.3'
        testGeneNames = strsplit(testGeneNames,split = ' ')[[1]]
        ifMarker = 0
        for(marker_i in 1:length(testGeneNames)){
            if(toupper(testGeneNames[marker_i]) %in% toupper(markerList)){
                ifMarker = ifMarker + 1
            }
        }
        if(ifMarker > 0 ){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}


check_subcellular_location = function(testLoc,locDescription){
    if(grepl(testLoc,locDescription)){
        return(1)
    }else{
        return(0)
    }
}

get_subcellular_location = function(ensembleID='ENSMUSG00000030735',
                                    fileFolder=paste0(adipocyte_dir,'MarkerPrediction/location/')){
    test = try(read.delim(file=paste0(fileFolder,ensembleID,'.txt'),
                          sep='\t',header = TRUE,stringsAsFactors = FALSE,check.names = FALSE))
    if(sum(grepl('Error',test))>0){return("")
    }else{
        tmpData = read.delim(file=paste0(fileFolder,ensembleID,'.txt'),
                             sep='\t',header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
        tmpLoc = c()
        if(sum(grepl('Mitochondrion',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Mitochondrion')
        }
        if(sum(grepl('Nucleus',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Nucleus')
        }
        if(sum(grepl('Cytoplasm',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Cytoplasm')
        }
        if(sum(grepl('Endoplasmic reticulum',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Endoplasmic reticulum')
        }
        if(sum(grepl('Lysosome',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Lysosome')
        }
        if(sum(grepl('Golgi',tmpData$`Subcellular location [CC]`))>0){
            tmpLoc=append(tmpLoc,'Golgi')
        }
        if(length(tmpLoc)>0){
            return(paste(tmpLoc,collapse = ';'))
        }else{
            return("")
        }
    }
}

## wrap text
wrap.it <- function(x, len){
    sapply(x, function(y) paste(strwrap(y, len), 
                                collapse = "\n"), 
           USE.NAMES = FALSE)
}
wrap.label <- function(x, len){
    if (is.list(x))
    {
        lapply(x, wrap.it, len)
    } else {
        wrap.it(x, len)
    }
}

plot_david_output = function(inFilePath=paste0(adipocyte_dir,'MarkerPrediction/david/bp.txt'),
                             figSavePath=NULL,splitSymbol='~',top=5,label.cex=0.6,...){
    tmpData = read.delim(file = inFilePath,sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE)
    if(nrow(tmpData)>top){
        tmpData = tmpData[1:top,]  
    }
    tmpData = tmpData[order(tmpData$PValue,decreasing = TRUE),]
    if(is.null(splitSymbol)){
        tmpData$TermsNumber = tmpData$Term
        tmpData$Terms =tmpData$Term
    }else{
        tmpData$TermsNumber = sapply(tmpData$Term,function(x) strsplit(x, split = splitSymbol)[[1]][1])
        tmpData$Terms = sapply(tmpData$Term,function(x) strsplit(x, split = splitSymbol)[[1]][2])
    }
    
    if(is.null(figSavePath)){
        tmpBar = barplot(-log(tmpData$PValue),horiz = TRUE,
                         #names.arg = wrap.label(tmpData$Terms,5),
                         beside = TRUE, las=1,
                         border = NA,
                         xlab = "-log(p-value)",
                         #xlim = c(0,round(max(-log(tmpData$PValue)),0)+1),
                         ...)
        text(x = -0.05, y = tmpBar, srt = 0,adj = 1, 
             labels = wrap.label(tmpData$Terms,15), 
             xpd = TRUE,cex = label.cex)
        #text(x = -log(tmpData$PValue)+0.2, y=tmpBar,label = tmpData$Count,col = 'black',cex = 0.5)
        
    }else{
        pdf(file=figSavePath)
        tmpBar = barplot(-log(tmpData$PValue),horiz = TRUE,
                         #names.arg = wrap.label(tmpData$Terms,5),
                         beside = TRUE, las=1,border = NA,
                         xlab = "-log(p-value)",
                         #xlim = c(0,round(max(-log(tmpData$PValue)),0)+1),
                         ...)
        text(x = -0.05, y = tmpBar, srt = 0,adj = 1, 
             labels = wrap.label(tmpData$Terms,15), 
             xpd = TRUE,cex = label.cex)
        #text(x = -log(tmpData$PValue)+0.2, y=tmpBar,label = tmpData$Count,col = 'black',cex = 0.5)
        #   # add label to x axis
        #   text(x = tmpBar, y = par("usr")[3] - 1, srt = 45,
        #        adj = 1, labels = tmpData$Terms, xpd = TRUE,cex = 0.35)
        
        dev.off()
    }
}

from_ensemblPeptide_to_genes = function(markersMapping=markerIDs,ensemblPeptideID,colName="ensembl_gene_id"){
    ensemblPeptideID = ifelse(grepl('\\.',ensemblPeptideID),strsplit(ensemblPeptideID,split = '\\.')[[1]][2],ensemblPeptideID)
    return(paste(unique(markersMapping[grepl(ensemblPeptideID,markersMapping$ensembl_peptide_id),colName]),collapse = ';'))
}

remove_redundant_terms = function(inFilePath=paste0(markerFolder,'david/batBPfatClassificationMedium.txt')){
    tmpLines =  read_files_by_line(file_path = inFilePath)
    tmpIndex = grep('Annotation Cluster',tmpLines)
    tmpIndex = append(tmpIndex,length(tmpLines))
    for(index_i in 1:(length(tmpIndex)-1)){
        if(index_i==1){
            tmpData = read.delim(file = inFilePath,sep = '\t',
                                 header = TRUE,stringsAsFactors = FALSE,
                                 skip = tmpIndex[index_i],nrows = (tmpIndex[1+index_i]-tmpIndex[index_i]-3))
            tmpData$cluster = index_i
            newData = tmpData[1,]
        }else{
            tmpData = read.delim(file = inFilePath,sep = '\t',
                                 header = TRUE,stringsAsFactors = FALSE,
                                 skip = tmpIndex[index_i],nrows = (tmpIndex[1+index_i]-tmpIndex[index_i]-3))
            tmpData$cluster = index_i
            newData = rbind(newData,tmpData[1,])
            
        }
    }
    write.table(newData,file = gsub('.txt$','Unique.txt',inFilePath),
                sep='\t',row.names = FALSE,quote = FALSE)
}


#### boxplot all genes using for
boxplotMarkers = function(markerData,mitoList=mitocarta,...){
  for(tmp_i in 1:(ncol(markerData)-1)){
    boxplot(markerData[,tmp_i]~markerData[,ncol(markerData)],
            data=markerData,main=colnames(markerData)[tmp_i],
            col=c('brown2','brown','papayawhip','orange'),
            col.main = ifelse(colnames(markerData)[tmp_i] %in% mitoList$Symbol,'blue4','deepskyblue2'),
            ylab='Gene expression (rlog)',...)
  }
}


##################################3333
###### functions for markerRegulatoryNetwork.R
#####################################
library(stringr)
split_string = function(a_string,split_symbol,aimed_index){
  # symbol is '|'
  # aimed_index is a number, which item you want to return after spliting
  temp_pattern = paste0('\\',split_symbol)
  id = str_split(a_string,pattern = temp_pattern)[[1]][aimed_index]
  return(id)
}

split_a_string_to_list = function(a_string,pattern){
  str_list = str_split(a_string,pattern = pattern)[[1]]
  str_list = str_list[str_list!= ""]
  return(str_list)
}

change1LineToMultipleLines = function(singleRowFromtfTable){
  regulatedGenes = split_a_string_to_list(singleRowFromtfTable$Target.genes,pattern = '\\,')
  tmpPPI = as.data.frame(matrix(data = NA,nrow = length(regulatedGenes),ncol = 1+ncol(singleRowFromtfTable)))
  colnames(tmpPPI) = c(colnames(singleRowFromtfTable),'regulatedGene')
  tmpPPI$X..Rank = singleRowFromtfTable$X..Rank
  tmpPPI$Motif.id = singleRowFromtfTable$Motif.id
  tmpPPI$AUC = singleRowFromtfTable$AUC
  tmpPPI$NES = singleRowFromtfTable$NES
  tmpPPI$ClusterCode = singleRowFromtfTable$ClusterCode
  tmpPPI$Transcription.factor = singleRowFromtfTable$Transcription.factor
  tmpPPI$Target.genes = singleRowFromtfTable$Target.genes
  tmpPPI$selectedTF = singleRowFromtfTable$selectedTF
  tmpPPI$regulatedGene = regulatedGenes
  return(tmpPPI)
}


getiRegulonPPI = function(iregulonFile = paste0(markerFolder,'iregulon/iRegulonTSS20kb.tsv'),
                          markerInformationFile,TFsList= c("Gata6","Ppara","Nr4a1","Ppargc1a")){
  iregulon = read.delim(file = iregulonFile,skip = 11,sep='\t',stringsAsFactors = FALSE)
  iregulon = iregulon[order(iregulon$NES,decreasing = TRUE),]
  clusterCodes = unique(iregulon$ClusterCode)
  tfTable = as.data.frame(matrix(data = NA,nrow = length(clusterCodes),ncol = 4))
  colnames(tfTable) = c('clusterCode','TFs','representativeTF','targetGenes')
  for(cluster_i in 1:length(clusterCodes)){
    tfTable$clusterCode[cluster_i] = clusterCodes[cluster_i]
    tmpTFs = unique(split_a_string_to_list(paste(iregulon$Transcription.factor[iregulon$ClusterCode==clusterCodes[cluster_i]],
                                                 collapse = ','),
                                           pattern = '\\,'))
    tfTable$TFs[cluster_i] = paste(tmpTFs,collapse = ',')
    tfTable$representativeTF[cluster_i] = tmpTFs[1]
    tmpTargets = unique(split_a_string_to_list(paste(iregulon$Target.genes[iregulon$ClusterCode==clusterCodes[cluster_i]],
                                                     collapse = ','),
                                               pattern = '\\,'))
    tfTable$targetGenes[cluster_i] = paste(tmpTargets,collapse = ',')
    rm(tmpTargets,tmpTFs)
  }
  tfTable = tfTable[!is.na(tfTable$representativeTF),] # remove NA
  tfTable = tfTable[!duplicated(tfTable$representativeTF),]
  write.table(tfTable,
              file = gsub('.tsv$','TFsClusters.txt',iregulonFile),
              sep = '\t',quote = FALSE,row.names = FALSE)
  iregulonPPI = as.data.frame(matrix(data = NA,0,ncol = 1+ncol(tfTable)))
  colnames(iregulonPPI)= c(colnames(tfTable),'regulatedGene')
  for(ppi_i in 1:nrow(tfTable)){
    iregulonPPI = rbind(iregulonPPI,change1LineToMultipleLines(tfTable[ppi_i,]))
  }
  write.table(iregulonPPI,
              file = gsub('.tsv$','PPI.txt',iregulonFile),
              sep = '\t',quote = FALSE,row.names = FALSE)
  ### ppi Attributes
  mouseMarker = read.delim(file = markerInformationFile,header = TRUE,stringsAsFactors = FALSE)
  ppiAttribute = as.data.frame(matrix(NA,
                                      nrow = length(unique(iregulonPPI$representativeTF)),
                                      ncol = 2))
  colnames(ppiAttribute) = c('GeneName','Type')
  ppiAttribute$GeneName = unique(iregulonPPI$representativeTF)
  ppiAttribute$Type = 'TF'
  ppiAttribute = rbind(ppiAttribute,mouseMarker[,c('GeneName','Type')])
  ppiAttribute$Class = ifelse(ppiAttribute$Type == 'TF','TF','Markers')
  write.table(ppiAttribute,
              file = gsub('.tsv$','PPIAttributes.txt',iregulonFile),
              sep = '\t',quote = FALSE,row.names = FALSE)

  return(iregulonPPI)
}



getiRegulonRegulatoryNetwork = function(iregulonFile = paste0(markerFolder,'iregulon/iRegulonTSS20kb.tsv'),
                          markerInformationFile,logFCcutoff=1.5,pvalueCutoff = 0.01,
                          sigFile=paste0(markerFolder,'df_exp.txt'),
                          mitoCartaList = mitocarta$Symbol){
  ## read iregulon
  
  iregulon = read.delim(file = iregulonFile,skip = 11,sep='\t',stringsAsFactors = FALSE)
  TFs = c()
  for(cluster_i in 1:nrow(iregulon)){
    tmpTFs = unique(split_a_string_to_list(iregulon$Transcription.factor[cluster_i],pattern = '\\,'))
    TFs = unique(append(TFs,tmpTFs))
  }
  
  
  ## read significant file
  sig = read.delim(file=sigFile,sep = '\t',stringsAsFactors = FALSE)
  sig$X = rownames(sig)
  sig = sig[(abs(sig$logFC)>logFCcutoff & sig$adj.P.Val<pvalueCutoff),]
  dim(sig)
  sig = merge(sig,getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                        filters = 'ensembl_gene_id',
                        values = sig$X,
                        mart = ensemblMartMouse),
              by.x='X',by.y='ensembl_gene_id')
  
  TFsList = sig$external_gene_name[(toupper(sig$external_gene_name) %in% toupper(TFs))]
  
  iregulon = iregulon[order(iregulon$NES,decreasing = TRUE),]
  iregulon$selectedTF = NA_character_
  for(m_i in 1:nrow(iregulon)){
    selectedtf = c()
    for(tf_i in 1:length(TFsList)){
      if(grepl(TFsList[tf_i],iregulon$Transcription.factor[m_i])){
        selectedtf = append(selectedtf,TFsList[tf_i])
      }
    }
    iregulon$selectedTF[m_i]=paste(selectedtf,collapse = ',')
  }
  iregulon = iregulon[iregulon$selectedTF !='',]
  iregulonNew = as.data.frame(matrix(NA,0,ncol(iregulon)))
  colnames(iregulonNew) = colnames(iregulon)
  ### split
  for(m_i in 1:nrow(iregulon)){
    if(grepl(',',iregulon$selectedTF[m_i])){
      tmpTFs = unique(split_a_string_to_list(a_string = iregulon$selectedTF[m_i],
                                      pattern = ','))
      tmpdf = as.data.frame(matrix(NA,nrow = length(tmpTFs),ncol(iregulon)))
      colnames(tmpdf) = colnames(iregulon)
      tmpdf[,1:7] = iregulon[m_i,1:7]
      tmpdf$selectedTF = tmpTFs
      iregulonNew = rbind(iregulonNew,tmpdf)
    }else{
      iregulonNew = rbind(iregulonNew,iregulon[m_i,])
    }
  }
  
  iregulonPPI = as.data.frame(matrix(data = NA,0,ncol = 1+ncol(iregulonNew)))
  colnames(iregulonPPI)= c(colnames(iregulonNew),'regulatedGene')
  for(ppi_i in 1:nrow(iregulonNew)){
    iregulonPPI = rbind(iregulonPPI,change1LineToMultipleLines(iregulonNew[ppi_i,]))
  }
  ###
  iregulonPPI$filter = apply(iregulonPPI[,c("selectedTF","regulatedGene")],1,function(x) paste(x,collapse = ','))
  iregulonPPI = iregulonPPI[!duplicated(iregulonPPI$filter),]
  iregulonPPI$line = iregulonPPI$selectedTF
  write.table(iregulonPPI,
              file = gsub('.tsv$',paste0('_logFC',logFCcutoff,'Pvalue',pvalueCutoff,'_PPI.txt'),iregulonFile),
              sep = '\t',quote = FALSE,row.names = FALSE)
  ### ppi Attributes
  mouseMarker = read.delim(file = markerInformationFile,header = TRUE,stringsAsFactors = FALSE)
  ppiAttribute = as.data.frame(matrix(NA,
                                      nrow = length(unique(iregulonPPI$selectedTF)),
                                      ncol = 3))
  colnames(ppiAttribute) = c('GeneName','Type','MitoCarta')
  ppiAttribute$GeneName = unique(iregulonPPI$selectedTF)
  ppiAttribute$Type = 'TF'
  ppiAttribute$MitoCarta = sapply(toupper(ppiAttribute$GeneName),function(x) ifelse(x %in% toupper(mitoCartaList),'Yes','No'))
  ppiAttribute = merge(ppiAttribute,sig[,c("logFC","external_gene_name")],all.x = TRUE,by.x = 'GeneName',by.y='external_gene_name')
  
  ppiAttribute = rbind(ppiAttribute,mouseMarker[,c('GeneName','Type','MitoCarta','logFC')])
  ppiAttribute$Class = ifelse(ppiAttribute$Type == 'TF','TF','Markers')
  write.table(ppiAttribute,
              file = gsub('.tsv$',paste0('_logFC',logFCcutoff,'Pvalue',pvalueCutoff,'_PPIAttributes.txt'),iregulonFile),
              sep = '\t',quote = FALSE,row.names = FALSE)
  return(iregulonPPI)
}

## for fig2c_pca_separate.R
add_0_to_study_number = function(studyNumber){
  number = str_split(studyNumber, pattern = 'M|R')[[1]][2]
  if(nchar(number) < 2){
    return(paste0('0',number))
  }else{
    return(number)
  }
}

plot_pca_separately = function(filePath, species = 'mouse'){
  tmpData =  read.delim(file = filePath,sep = '\t',header = TRUE,
                        row.names = 1,check.names = FALSE)
  tmpDataPCA=prcomp(t(tmpData))
  screeplot(tmpDataPCA)
  pcaDat = as.data.frame(tmpDataPCA$x[,1:3])
  tmpDataPCAsummary = summary(tmpDataPCA)
  # add information to pca result
  pcaDat$study = as.factor(str_extract(rownames(pcaDat), pattern = '[M|R][0-9]+'))
  pcaDat$study = sapply(pcaDat$study, function(x) paste0(ifelse(grepl('M',x),'M','R'), add_0_to_study_number(x)))
  if(species == 'human'){
    pcaDat$study = paste0('h',pcaDat$study)
  }
  pcaDat$tissue = as.factor(unname(sapply(rownames(pcaDat),function(x) get_adipocyte_type(x))))
  #unique(get_adipocyte_colors(pcaDat$tissue))
  if(sum(grepl('BAT\\(T\\)', pcaDat$tissue)) >0 ){
    pca_colors = c("brown2", "brown", "papayawhip", "orange")
  }else {
    pca_colors = c("brown2", "papayawhip", "orange")
  }
  
  # scatter plot
  ggplot(pcaDat, aes(x = PC1, y =  PC2, color = tissue)) + 
    geom_point(size = 1.5) +
    facet_grid(study ~ .) +
    scale_color_manual(values = pca_colors ) + 
    theme(panel.background = element_rect(fill = 'lightgray'), # set background color
          panel.grid = element_blank(), # remove the background grid
          legend.key = element_rect(fill = 'lightgray'), # change background color of legend
          strip.background = element_rect(fill = 'white'), # change the background color of each row 'M01'
          axis.text.y = element_text(size = 5)) + # change y label font size
    labs(x = paste0('PC1 (',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC1'],2)*100,'%)'),
         y = paste0('PC2 (',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC2'],2)*100,'%)'),
         color = "")
  
}

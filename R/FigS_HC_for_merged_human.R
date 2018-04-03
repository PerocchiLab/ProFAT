# replot the HC for microArray and RNA-Seq data
# jiang
# 2017/09/25

library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(car)
library(ggplot2)
library(randomcoloR)
library(xlsx)
library(stringi)
library(stringr)
library(hashmap)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))


samples = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/human_datasets.xls'), 
                    sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
samples = samples[,1:5]

####################
### map latest sample ID to old sample IDs
old_sample_id = read.delim(file = paste0(adipocyte_dir, 'DataInfo/NoCore_Dataset.txt'),
                           header = TRUE, stringsAsFactors = FALSE)
old_sample_id = old_sample_id[grepl('^h', old_sample_id$DataSetID), c(1,2,6)]

sample_id = merge(old_sample_id, samples, by.x = 'SampleID',by.y = 'Sample.ID')
sample_id_map = hashmap(sample_id$StudyID, sample_id$Dataset.ID)


geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.human.All.txt'),
                        paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.human.All.rlog.txt'))
#############################
## HC microarray
dist_metric = "euclidean"
geneExpressionFile = geneExpressionFiles[1]
saveFile = gsub("txt$", paste0("pheatmap.euclidean.pdf"), geneExpressionFile)

tmpData =  read.delim(file=geneExpressionFile,
                      sep = '\t',header = TRUE,row.names = 1,check.names = FALSE)

colors = colorRampPalette( rev(brewer.pal(8, "Greens")) )(255)
font_size = ceiling(300/dim(tmpData)[2])
if (font_size>10){
  font_size = 10
}

#datDist  = dist(t(tmpData),method = "euclidean")
datDist = Dist(t(tmpData), method=dist_metric)
datDistMatrix = as.matrix(datDist)
# rowname original format 'BAT.line2_hR2'
rownames(datDistMatrix)= sample_id_map[[rownames(datDistMatrix)]]
colnames(datDistMatrix)=NULL
  pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
  pheatmap(datDistMatrix,
           clustering_distance_rows=datDist,
           clustering_distance_cols=datDist,
           fontsize_col = font_size,
           col=colors, fontsize_row=font_size,
           main='euclidean', border_color	=NA)
  dev.off()
  
  
  ##############################333
  ## HC for RNA-Seq
  dist_metric = "euclidean"
  geneExpressionFile = geneExpressionFiles[2]
  saveFile = gsub("txt$", paste0("pheatmap.euclidean.pdf"), geneExpressionFile)
  
  tmpData =  read.delim(file=geneExpressionFile,
                        sep = '\t',header = TRUE,row.names = 1,check.names = FALSE)
  
  colors = colorRampPalette( rev(brewer.pal(8, "Greens")) )(255)
  font_size = ceiling(300/dim(tmpData)[2])
  if (font_size>10){
    font_size = 10
  }
  
  #datDist  = dist(t(tmpData),method = "euclidean")
  datDist = Dist(t(tmpData), method=dist_metric)
  datDistMatrix = as.matrix(datDist)
  # rowname original format 'BAT.line2_hR2'
  rownames(datDistMatrix)= sample_id_map[[rownames(datDistMatrix)]]
  colnames(datDistMatrix)=NULL
  pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
  pheatmap(datDistMatrix,
           clustering_distance_rows=datDist,
           clustering_distance_cols=datDist,
           fontsize_col = font_size,
           col=colors, fontsize_row=font_size,
           main='euclidean', border_color	=NA)
  dev.off()

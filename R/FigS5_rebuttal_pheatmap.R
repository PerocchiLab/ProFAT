# caculate and display the distance among samples within 1 study
# jiang
# 2017/01/26


# ## copy this script from jiang's dir to share/dir
# system('cp /data/home/jiang/projects/scripts/adipocyte/20170126/* /data/home/share/Projects/Adipocyte/rCode/')

library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(car)
library(foreach)
library(doParallel)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
# code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))



dataFile = paste0(adipocyte_dir,'MicroArray/Old_All_BeBrW/rbatch.All.ND.13.txt')
figurePath = gsub("txt$", paste0("pheatmap.euclidean.sidecolorbar.pdf"), dataFile)
geneExpressionFile = dataFile
saveFile = figurePath
dist_metric = "euclidean"

tmpData =  read.delim(file=geneExpressionFile,
                      sep = '\t',header = TRUE,row.names = 1,check.names = FALSE)

colors = colorRampPalette( rev(brewer.pal(8, "Greens")) )(255)
font_size = ceiling(300/dim(tmpData)[2])
if (font_size>10){
  font_size = 10
}


datDist = Dist(t(tmpData), method=dist_metric)
datDistMatrix = as.matrix(datDist)
# rowname original format 'BAT.line2_hR2'
rownames(datDistMatrix)= sapply(rownames(datDistMatrix),
                                function(x) ifelse(grepl('_hM|_hR',x),gsub('\\.','\\_',x),x))

colnames(datDistMatrix) = rownames(datDistMatrix)
annotation_row = data.frame(GeneClass = 
                              factor(ifelse(grepl('iWAT\\(T\\)\\_M11\\(r2\\)', rownames(datDistMatrix)), 'remove','keep' )))
rownames(annotation_row) = rownames(datDistMatrix)
ann_colors = list(GeneClass = c(remove = "red", keep = "black"))

pdf(file = saveFile, onefile = FALSE, width = (18.3/2.54), height = (24.7/2.54))
pheatmap(datDistMatrix,
         clustering_distance_rows=datDist,
         clustering_distance_cols=datDist,
         fontsize_col = font_size,
         col=colors, fontsize_row=font_size,
         show_rownames = TRUE, show_colnames = FALSE,
         annotation_col = annotation_row, annotation_row = annotation_row, 
         annotation_colors = ann_colors, annotation_legend = FALSE)
dev.off()
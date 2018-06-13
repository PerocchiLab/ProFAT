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
###################################################
## pheatmap
###################################################
 cl = makeCluster(30)
 registerDoParallel(cl)


p_forms<-c("MicroArray", "RNASeq")
#p_forms<-"MicroArray"
for (t_platform in p_forms){
  arrayExpressFiles = list.files(path = paste0(adipocyte_dir,t_platform),
                                 pattern = "(ND|rlog).*.txt",recursive = TRUE,
                                 full.names = TRUE)
  arrayExpressFiles = arrayExpressFiles[!grepl('old',arrayExpressFiles)]

  #dist_methods<-c("euclidean", "pearson", "spearman")
 # dist_methods<-c("euclidean","showGeneValue")
  #dist_methods<-c("euclidean","showGeneValue","showGeneValueScaleColumn")
  dist_methods<-c("euclidean")
  
  
  foreach(i=1:length(arrayExpressFiles),.packages = c('pheatmap','amap','RColorBrewer'),.export = c('plots_heatmap_distance')) %dopar% {
   # for(i in 1:length(arrayExpressFiles)){
    #   tmpPlotsDir = create_new_dir(create_new_dir(dirname(arrayExpressFiles[i])),'/plots/')
    #   tmpSaveFile = paste0(tmpPlotsDir,basename(dirname(arrayExpressFiles[i])),'_heatmap_distance.pdf')
    for (tmp_metric in dist_methods){
      tmpSaveFile = gsub("txt$", paste0("pheatmap.", tmp_metric, ".pdf"), arrayExpressFiles[i])
      tryPlot = try(
        plots_heatmap_distance(geneExpressionFile = arrayExpressFiles[i],
                               saveFile = tmpSaveFile, dist_metric = tmp_metric, main=tmp_metric,
                               border_color	=NA)
      )
      if(class(tryPlot) != 'try-error'){
        plots_heatmap_distance(geneExpressionFile = arrayExpressFiles[i],
                               saveFile = tmpSaveFile, dist_metric = tmp_metric, main=tmp_metric,
                               border_color	=NA)
      }

    }
    #       plots_pca(geneExpressionFile = arrayExpressFiles[i],
    #                 saveFile = gsub('txt$','pca.pdf',arrayExpressFiles[i]),
    #                 colorGroups = NULL,
    #                 pchGroups=NULL,
    #                 main=paste0(gsub(".txt$", "", basename(arrayExpressFiles[i])), " PCA"))
  }
}


stopCluster(cl)

# tmp_metric = "showGeneValueScaleColumn"
# rnaFile = 
# tmpSaveFile = gsub("txt$", paste0("pheatmap.", tmp_metric, ".pdf"), rnaFile)
# plots_heatmap_distance(geneExpressionFile = rnaFile,
#                        saveFile = tmpSaveFile, 
#                        dist_metric = tmp_metric,
#                        main=tmp_metric,
#                        border_color	=NA)
# arrrayFile = 
# tmpSaveFile = gsub("txt$", paste0("pheatmap.", tmp_metric, ".pdf"), arrrayFile)
# plots_heatmap_distance(geneExpressionFile = arrrayFile,
#                        saveFile = tmpSaveFile, 
#                        dist_metric = tmp_metric,
#                        main=tmp_metric,
#                        border_color	=NA)
  
#################################################################################
## PCA
##########################################################################
####################################################################
## pca via plot_pca function
######################
# ## microarray
# dataFile = paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.txt')
# figurePath = gsub("txt$", "pca.pdf", dataFile)
# plots_pca(geneExpressionFile = dataFile,
#           saveFile = figurePath,
#           colorGroups = NULL,pchGroups = NULL,
#           ypc = 2,cex=1,radius=2)
# 
# ## RNAseq
# dataFile = paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.rlog.txt')
# figurePath = gsub("txt$", "pca.pdf", dataFile)
# plots_pca(geneExpressionFile = dataFile,
#           saveFile = figurePath,
#           colorGroups = NULL,pchGroups = NULL,
#           ypc = 2,cex=1,radius=1.7)

############################################################################
#dataFile = paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.agg.txt')
dataFile = paste0(adipocyte_dir,'MicroArray/Old_All_BeBrW/rbatch.All.ND.13.txt')

figurePath = gsub("txt$", paste0("pheatmap.euclidean.pdf"), dataFile)
plots_heatmap_distance(geneExpressionFile = dataFile,
                       saveFile = figurePath,dist_metric = "euclidean", 
                       main='euclidean', border_color	=NA)


### hierachical clustering for all human dataset
geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.human.All.txt'),
                        paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.human.All.rlog.txt'))

figurePath = gsub("txt$", paste0("pheatmap.euclidean.pdf"), geneExpressionFiles[1])
plots_heatmap_distance(geneExpressionFile = geneExpressionFiles[1],
                       saveFile = figurePath,dist_metric = "euclidean", 
                       main='euclidean', border_color	=NA)

figurePath = gsub("txt$", paste0("pheatmap.euclidean.pdf"), geneExpressionFiles[2])
plots_heatmap_distance(geneExpressionFile = geneExpressionFiles[2],
                       saveFile = figurePath,dist_metric = "euclidean", 
                       main='euclidean', border_color	=NA)

#
# dataFile = paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.agg.txt')
# figurePath = gsub("txt$", paste0("pheatmap.euclidean.greys.pdf"), dataFile)
# plots_heatmap_distance(geneExpressionFile = dataFile,
#                        saveFile = figurePath,dist_metric = "euclidean", main='euclidean')

## change the colors in function plots_heatmap_distance
# dataFile = paste0(adipocyte_dir,'MicroArray/old_All_BeBrW/rbatch.All.ND.txt')
# figurePath = gsub("txt$", paste0("pheatmap.euclidean.ColNames.pdf"), dataFile)
# plots_heatmap_distance(geneExpressionFile = dataFile,
#                        saveFile = figurePath,
#                        dist_metric = "euclidean", main='euclidean')


# ## distance among RNASeq samples in each study
# seqExpressFiles = list.files(path = paste0(adipocyte_dir,'RNASeq/'),
#                                pattern = 'fpkm.txt',recursive = TRUE,
#                                full.names = TRUE)
# for(i in 1:length(seqExpressFiles)){
#   #tmpPlotsDir = create_new_dir(create_new_dir(dirname(seqExpressFiles[i])),'/plots/')
#   tmpSaveFile = gsub("txt$", "pheatmap.pdf", seqExpressFiles[i])
#   plots_heatmap_distance(geneExpressionFile = seqExpressFiles[i],
#                          saveFile = tmpSaveFile)
# }
#
#
# ## RNAseq
# dataFile = paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rlog.seqcount.txt')
# figurePath = gsub("txt$", "pheatmap.pdf", dataFile)
# plots_heatmap_distance(geneExpressionFile = dataFile,
#                        saveFile = figurePath,
#                        fontsize_row=6)
#
# dataFile = paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/nc.seqcount.txt')
# figurePath = gsub("txt$", "pheatmap.pdf", dataFile)
# plots_heatmap_distance(geneExpressionFile = dataFile,
#                        saveFile = figurePath,
#                        fontsize_row=6)

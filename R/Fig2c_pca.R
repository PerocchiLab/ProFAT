# caculate and display the distance among samples within 1 study
# jiang
# 2017/02/06
# ## copy this script from jiang's dir to share/dir
# system('cp /data/home/jiang/projects/scripts/adipocyte/20170126/* /data/home/share/Projects/Adipocyte/rCode/')

library("pheatmap")
library("RColorBrewer")
library("amap")
library('rafalib')
library("car")

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))
#####################################################################3
## test PCA
##################################

# geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.txt'),
#                         paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.rlog.txt'))

geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.13.txt'))

for(i in 1:length(geneExpressionFiles)){
  if(grepl('MicroArray',geneExpressionFiles[i])){
    preDefinePch = rev(c(0:10,12:13))
  }else{
    preDefinePch = c(0:10,12:13)
  }


  for(i in seq(1.5,1.5,0.5)){
    #saveFile = paste0(plots_dir,basename(saveFile))

    colorGroups = NULL;pchGroups=NULL;ypc = 2;
    #### test pca
    tmpData =  read.delim(file=geneExpressionFiles[i],sep = '\t',header = TRUE,
                          row.names = 1,check.names = FALSE)
    # if the sample group information is provided by the author,we use the given informaton,
    # otherwise we extract the information from the colnames of the given dataset
    if(is.null(colorGroups)){
      tmpColorGroups = unname(sapply(colnames(tmpData),function(x) get_adipocyte_type(x)))
    }else{
      tmpColorGroups = colorGroups
    }
    if(is.null(pchGroups)){
      tmpPchGroups = str_extract(colnames(tmpData), pattern = 'M[0-9]+')
    }else{
      tmpPchGroups = pchGroups
    }

    tmpDataPCA=prcomp(t(tmpData))
    tmpDataPCAsummary = summary(tmpDataPCA)
    ###add sample type to tmpDataPCA$x
    dfCircle = as.data.frame(tmpDataPCA$x)
    dfCircle$type = as.factor(sapply(row.names(tmpDataPCA$x),function(x) get_adipocyte_type(x)))
    #pdf(file = saveFile,onefile = FALSE)
    setEPS()
    #saveFile = gsub("txt$", paste0("pca_size",i,".pdf"), geneExpressionFiles[i])
    saveFile = gsub("txt$", paste0("pca_size",i,".eps"), geneExpressionFiles[i])
    postscript(saveFile)
    #par(bg='gray')
    plot(0, 0, type="n", ann=FALSE, axes=FALSE)
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col="gray", border=NA)
    par(new=TRUE)
    plot(tmpDataPCA$x[,1],tmpDataPCA$x[,2],
         xlim = c(min(tmpDataPCA$x[,1]-2),max(tmpDataPCA$x[,1]+15)),
         col=get_adipocyte_colors(tmpColorGroups),
         pch=preDefinePch[as.fumeric(tmpPchGroups)],cex=i,
         xlab = paste0('PC1: ',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC1'],2)*100,'% variance'),
         ylab = paste0('PC',ypc,': ',round(tmpDataPCAsummary$importance["Proportion of Variance",][paste0('PC',ypc)],2)*100,'% variance'),
         add=TRUE)


    #   # add legend by groups, on the topleft of the figure
    legend('bottom',unique(tmpColorGroups),col = get_adipocyte_colors(unique(tmpColorGroups)),
           lty = 1,bty = 'n',text.col = get_adipocyte_colors(unique(tmpColorGroups)))
    legend('right',unique(tmpPchGroups),col = 'black',bty = 'n',
           pch = unique(preDefinePch[as.fumeric(tmpPchGroups)]),
           cex = 1)
    dev.off()
  }

}


##########################################333333333
###v  PCA on Martin's data
###############################################
inhouseMarkerExpressionFile = paste0(adipocyte_dir,'RNASeq/m00R_BrBeW/tg/rlog.txt')
plots_pca(geneExpressionFile = inhouseMarkerExpressionFile,
          saveFile = gsub('txt$','pca.pdf',inhouseMarkerExpressionFile))
plots_pca(geneExpressionFile = inhouseMarkerExpressionFile,
          saveFile = paste0(plots_dir,'test.inhouse.pca.pdf'))

## PCA for all 13 microArray studies
microArray = paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.13.txt')
plots_pca(geneExpressionFile = microArray,
          saveFile = gsub('txt$','pca.pdf',geneExpressionFiles), cex = 1.5, lwd=1.5)
########################################################################

microArrayFile = list.files(path = paste0(adipocyte_dir, 'MicroArray/'), 
                            pattern = 'ND_[h|m][0-9]+M.*.txt', 
                            recursive = TRUE, full.names = TRUE)
for(single_array in microArrayFile){
  plots_pca(geneExpressionFile = single_array,
            saveFile = gsub('txt$','pca.pdf',single_array),
            cex = 1.5, lwd = 2)
}

seqFile = list.files(path = paste0(adipocyte_dir, 'RNASeq/'), 
                            pattern = 'rlog.txt', 
                            recursive = TRUE, full.names = TRUE)
seqFile = seqFile[!grepl('old_', seqFile)]
for(single_seq in seqFile){
  plots_pca(geneExpressionFile = single_seq,
            saveFile = gsub('txt$','pca.pdf',single_seq),
            cex = 1.5, lwd = 2)
}

#######################################################################

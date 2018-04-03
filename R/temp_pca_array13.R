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
#################################################################################

get_adipocyte_type = function(sampleName){
  if(grepl('BAT\\(T\\)',sampleName)){
    return('BAT(T)')
  }else if(grepl('BAT',sampleName)){
    return('BAT')
  }else if(grepl('WAT\\(T\\)',sampleName)){
    return('WAT(T)')
  }else if(grepl('WAT\\(T\\_', sampleName)){
    return(str_split(sampleName, pattern = '_M[0-9]+')[[1]][1])
  }else {
    return('WAT')
  }
}


get_adipocyte_colors = function(adipocyte_character){
  adipocyte_character = gsub('BAT\\(T\\)','brown',adipocyte_character)
  adipocyte_character = gsub('BAT$','brown2',adipocyte_character)
  adipocyte_character = gsub('.*WAT\\(T\\)','orange',adipocyte_character)
  adipocyte_character = gsub('.*WAT$','papayawhip',adipocyte_character)
  adipocyte_character = gsub('iWAT\\(T\\_CL)','green',adipocyte_character)
  adipocyte_character = gsub('iWAT\\(T\\_RS)','yellow',adipocyte_character)
  adipocyte_character = gsub('iWAT\\(T\\_RG)','blue',adipocyte_character)
  return(adipocyte_character)
}
######################################################################################


microArrayFile = list.files(path = paste0(adipocyte_dir, 'MicroArray/'), 
                            pattern = 'ND_[h|m][0-9]+M.*.txt', 
                            recursive = TRUE, full.names = TRUE)
array13 = microArrayFile[grepl('m13M', microArrayFile)] # the 13th microArray only

plots_pca(geneExpressionFile = array13,
          saveFile = gsub('txt$','13pca.pdf',array13), cex = 1.5, lwd=1.5)


## PCA for all 13 microArray studies
microArray = paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.13.txt')
plots_pca(geneExpressionFile = microArray,
          saveFile = gsub('txt$','pca_distinct_Farmer.pdf',microArray), cex = 1.5, lwd=1.5)


preDefinePch = c(0:10,12:13)
colorGroups = NULL;pchGroups=NULL;ypc = 2;
#### test pca
tmpData =  read.delim(file=geneExpressionFiles,sep = '\t',header = TRUE,
                      row.names = 1,check.names = FALSE)
# if the sample group information is provided by the author,we use the given informaton,
# otherwise we extract the information from the colnames of the given dataset
if(is.null(colorGroups)){
  tmpColorGroups = unname(unlist(sapply(colnames(tmpData),function(x) get_adipocyte_type(x))))
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
saveFile = gsub("txt$",".eps", array13)
postscript(saveFile)
#par(bg='gray')
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
u <- par("usr") # The coordinates of the plot area
rect(u[1], u[3], u[2], u[4], col="gray", border=NA)
par(new=TRUE)
plot(tmpDataPCA$x[,1],tmpDataPCA$x[,2],
     xlim = c(min(tmpDataPCA$x[,1]-2),max(tmpDataPCA$x[,1]+15)),
     col=get_adipocyte_colors(tmpColorGroups),
     pch=preDefinePch[as.fumeric(tmpPchGroups)],cex=1.5, lwd = 2,
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

# replot the PCA for microArray and RNA-Seq data
# jiang
# 2017/07/24

library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(car)
library(ggplot2)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
#code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))
#####################################################################3
## test PCA
##################################
# 
geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.All.ND.13.txt'),
                        paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.rlog.txt'))

# geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.human.All.txt'),
#                         paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.human.All.rlog.txt'))
### plot pca 
### plot PCA on microArray gene expression
plot_pca_separately(filePath = geneExpressionFiles[1])
ggsave(gsub('txt$', 'pca_separately.pdf', geneExpressionFiles[1]), 
       width = 15,height = 12,units = 'cm')

### plot PCA on RNAseq gene expression
plot_pca_separately(filePath = geneExpressionFiles[2])
ggsave(gsub('txt$', 'pca_separately.pdf', geneExpressionFiles[2]),
       width = 15,height = 6,units = 'cm')


##################################################
### show different colors on treated WAT from M13
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

tmpData =  read.delim(file = geneExpressionFiles[1],sep = '\t',header = TRUE,
                      row.names = 1,check.names = FALSE)
tmpDataPCA=prcomp(t(tmpData))
screeplot(tmpDataPCA)
dim(tmpDataPCA$x)
pcaDat = as.data.frame(tmpDataPCA$x[,1:3])
tmpDataPCAsummary = summary(tmpDataPCA)
# add information to pca result
pcaDat$study = as.factor(str_extract(rownames(pcaDat), pattern = '[M|R][0-9]+'))
pcaDat$study = sapply(pcaDat$study, function(x) paste0(ifelse(grepl('M',x),'M','R'), add_0_to_study_number(x)))
pcaDat$tissue = as.factor(unname(sapply(rownames(pcaDat),function(x) get_adipocyte_type(x))))
#unique(get_adipocyte_colors(pcaDat$tissue))

# scatter plot
ggplot(pcaDat, aes(x = PC1, y =  PC2, color = tissue)) + 
  geom_point(size = 1.5) +
  facet_grid(study ~ .) +
  scale_color_manual(values = c("brown2", "brown", "green", "blue", "yellow", "papayawhip", "orange")) + 
  theme(panel.background = element_rect(fill = 'lightgray'), # set background color
        panel.grid = element_blank(), # remove the background grid
        legend.key = element_rect(fill = 'lightgray'), # change background color of legend
        strip.background = element_rect(fill = 'white'), # change the background color of each row 'M01'
        axis.text.y = element_text(size = 5)) + # change y label font size
  labs(x = paste0('PC1 (',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC1'],2)*100,'%)'),
       y = paste0('PC2 (',round(tmpDataPCAsummary$importance["Proportion of Variance",]['PC2'],2)*100,'%)'),
       color = "")
ggsave(gsub('txt$', 'pca_separately_distinguish_M13.pdf', geneExpressionFiles[1]),
       width = 15,height = 15,units = 'cm')


# adipocyte_character = gsub('iWAT\\(T\\_CL)','green',adipocyte_character)
# adipocyte_character = gsub('iWAT\\(T\\_RS)','yellow',adipocyte_character)
# adipocyte_character = gsub('iWAT\\(T\\_RG)','blue',adipocyte_character)





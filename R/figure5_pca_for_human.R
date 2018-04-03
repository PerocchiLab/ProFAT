# replot the PCA for microArray and RNA-Seq data
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
library(oligo)
library(affxparser)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))


samples = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/human_datasets.xls'), 
                    sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
samples = samples[,1:5]
stri_extract_all(str = samples$Dataset.ID[92], 
                 regex = '[S|s]ubcutaneous neck|[D|d]eep neck|([W|B]AT.+forskolin)|sWAT[^\\_]*|PSC[^\\_]+|([B|W]AT[ preadipocytes]*)')[[1]]
samples$sample = sapply(samples$Dataset.ID, function(x) stri_extract_all(str = x, 
                                                                         regex = '[S|s]ubcutaneous neck|[D|d]eep neck|([W|B]AT.+forskolin)|sWAT[^\\_]*|PSC[^\\_]+|([B|W]AT[ preadipocytes]*)')[[1]])
samples$sample[grepl('subcutaneous neck', samples$sample)] = 'Subcutaneous neck'
samples$sample[grepl('deep neck', samples$sample)] = 'Deep neck'

table(samples$sample)
length(unique(samples$sample))

potential_colors = c('blue', 'darkslategray2','lightcoral','plum2','chartreuse3',
                     'darkgreen', 'deeppink1', 'gold3', 'yellow','violetred4', 
                     'red', 'cyan4', 'darkorange','midnightblue')

tissue_color = hashmap(unique(samples$sample), potential_colors)

####################
### map latest sample ID to old sample IDs
old_sample_id = read.delim(file = paste0(adipocyte_dir, 'DataInfo/NoCore_Dataset.txt'),
                           header = TRUE, stringsAsFactors = FALSE)
old_sample_id = old_sample_id[grepl('^h', old_sample_id$DataSetID), c(1,2,6)]

sample_id = merge(old_sample_id, samples, by.x = 'SampleID',by.y = 'Sample.ID')
sample_id_map = hashmap(sample_id$StudyID, sample_id$sample)


geneExpressionFiles = c(paste0(adipocyte_dir,'MicroArray/All_BeBrW/rbatch.human.All.txt'),
                        paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/rbatch.human.All.rlog.txt'))


#########################################################3
## microarray
###################
filePath = geneExpressionFiles[1]
species = 'human'
platform = 'array'

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
# pcaDat$tissue = as.factor(unname(sapply(rownames(pcaDat),function(x) get_adipocyte_type(x))))
pcaDat$tissue = sample_id_map[[rownames(pcaDat)]]
# #unique(get_adipocyte_colors(pcaDat$tissue))
# if(sum(grepl('BAT\\(T\\)', pcaDat$tissue)) >0 ){
#   pca_colors = c("brown2", "brown", "papayawhip", "orange")
# }else {
#   pca_colors = c("brown2", "papayawhip", "orange")
# }

pca_colors = unique(tissue_color[[pcaDat$tissue]])

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
ggsave(gsub('txt$', 'pca_separately.pdf', geneExpressionFiles[1]),
       width = 15,height = 6,units = 'cm')

########################3
## RNA-Seq
filePath = geneExpressionFiles[2]
species = 'human'

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
# pcaDat$tissue = as.factor(unname(sapply(rownames(pcaDat),function(x) get_adipocyte_type(x))))
pcaDat$tissue = sample_id_map[[rownames(pcaDat)]]
# #unique(get_adipocyte_colors(pcaDat$tissue))
# if(sum(grepl('BAT\\(T\\)', pcaDat$tissue)) >0 ){
#   pca_colors = c("brown2", "brown", "papayawhip", "orange")
# }else {
#   pca_colors = c("brown2", "papayawhip", "orange")
# }

pca_colors = unique(tissue_color[[pcaDat$tissue]])

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
ggsave(gsub('txt$', 'pca_separately.pdf', geneExpressionFiles[2]),
       width = 15,height = 8,units = 'cm')

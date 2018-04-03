# caculate and display the distance among samples within 1 study
# jiang
# 2017/07/26
# ## copy this script from jiang's dir to share/dir

library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(sva)
library(gplots)
library(xlsx)


adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))


## female mice
femaleData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m00R_BrBeW/tg/rlog.txt'),
                        sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
femaleData = femaleData[,grepl('\\(T\\)',colnames(femaleData))] # remove treated samples
femaleData$Ensembl = rownames(femaleData)

# male mouse
maleData = read.delim(file = paste0(adipocyte_dir,'RNASeq/m20R_BrW/tg/rlog.txt'),
                      sep = '\t',stringsAsFactors = FALSE,check.names = FALSE)
colnames(maleData) = gsub('_', '(T)_', colnames(maleData))
maleData$Ensembl = rownames(maleData)

## merge
mergedDat = merge(maleData, femaleData, by = 'Ensembl', all = FALSE)
rownames(mergedDat) = mergedDat$Ensembl
mergedDat = mergedDat[, !(grepl('Ensembl',colnames(mergedDat)))]
colnames(mergedDat) = gsub('_R0', '_R00', colnames(mergedDat))
batch_VM = ifelse(grepl('R20',colnames(mergedDat)), 'R20', 'R00')
mergedDat = ComBat(dat = mergedDat, batch = batch_VM)

#################################################################################################
### heatmap for markers expression in both male and female
mouseMarker = read.delim(file = paste0(adipocyte_dir,'BrW_Markers/marker.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Ensembl = rownames(mouseMarker)
mouseMarkerBat = read.delim(file = paste0(adipocyte_dir,'BrW_Markers/marker.BAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarkerWat = read.delim(file = paste0(adipocyte_dir,'BrW_Markers/marker.WAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Type = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'BAT','WAT'))


mergedDatMarker = mergedDat[rownames(mergedDat) %in% rownames(mouseMarker),]
mergedDatMarker$Ensembl = rownames(mergedDatMarker)
mergedDatMarker = merge(mouseMarker[,c('Ensembl',"logFC","Type", 'GeneName')], mergedDatMarker,
                        by = 'Ensembl', all = FALSE)
rownames(mergedDatMarker) = mergedDatMarker$GeneName
mergedDatMarker = mergedDatMarker[order(mergedDatMarker$Type,-(abs(mergedDatMarker$logFC))),]

mergedDatMarker = mergedDatMarker[,-c(1:4)]
#######################################
## plot heatmap with original gene expression

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markers_female_male_rlog.pdf'),
    onefile = FALSE,family = "Helvetica",width = 11,height = 10)
heatmap.2(as.matrix(mergedDatMarker),
          main = "Markers Expression", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          #margins = c(12,9),     # widens margins around plot
          # col = colorRampPalette(c('lightyellow2','blue'))(n = 50),       # use on color palette defined earlier
          col = colorRampPalette(c('darkgreen','yellow','red'))(n = 60),       # use on color palette defined earlier
          #col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
          #col = colorRampPalette(c("green", "orange", "firebrick3"))(50),
          #labRow = TRUE,
          keysize = 1.5, srtCol = 45,
          cexRow = 0.8, cexCol = 0.8, scale = 'none',
          #Rowv = as.dendrogram(row_cluster), # apply default clustering method
          #Colv = as.dendrogram(col_cluster),
          Rowv=NULL, dendrogram = 'col')
legend("topright",
       legend=c('R20: Male','R00: Female'),
       fill=c('white'),
       border=FALSE, bty="n", y.intersp = 1, cex=1.5)
dev.off()

## plot heatmap with original gene expression rlog (scale by row)

pdf(file = paste0(adipocyte_dir,'BrW_Markers/markers_female_male_rlog_scale.pdf'),
    onefile = FALSE,family = "Helvetica",width = 11,height = 10)
heatmap.2(as.matrix(mergedDatMarker),
          main = "Markers Expression", # heat map title
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          #margins = c(12,9),     # widens margins around plot
          # col = colorRampPalette(c('lightyellow2','blue'))(n = 50),       # use on color palette defined earlier
          col = colorRampPalette(c('darkgreen','yellow','red'))(n = 60),       # use on color palette defined earlier
          #col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
          #col = colorRampPalette(c("green", "orange", "firebrick3"))(50),
          #labRow = TRUE,
          keysize = 1.5, srtCol = 45,
          cexRow = 0.8, cexCol = 0.8, scale = 'row',
          #Rowv = as.dendrogram(row_cluster), # apply default clustering method
          #Colv = as.dendrogram(col_cluster),
          Rowv=NULL, dendrogram = 'col')
legend("topright",
       legend=c('R20: Male','R00: Female'),
       fill=c('white'),
       border=FALSE, bty="n", y.intersp = 1, cex=1.5)
dev.off()

###########################################################################
# two-way anova
dat = mergedDatMarker
tissue  = as.factor(unname(sapply(colnames(mergedDatMarker), function(x) ifelse(grepl('BAT',x), "BAT(T)", 'WAT(T)'))))
gender  = as.factor(unname(sapply(colnames(mergedDatMarker), function(x) ifelse(grepl('R20',x), "Male", 'Female'))))
dat = t(dat)
fit = manova(dat ~ tissue + gender)
summaryRes = summary.aov(fit)
summaryRes[[1]]
rownames(summaryRes[[1]])
summaryRes[[1]]$`Pr(>F)`[2]

## add pvalue to mergedDatMarker
mergedDatMarker = round(mergedDatMarker, digits = 4)
mergedDatMarker$`pvalue_for_gender` = 1
mergedDatMarker$`pvalue_for_tissue_type` = 1
for(i_gene in 1:nrow(mergedDatMarker)){
  mergedDatMarker$`pvalue_for_tissue_type`[i_gene] = summaryRes[[i_gene]]$`Pr(>F)`[1]
  mergedDatMarker$`pvalue_for_gender`[i_gene] = summaryRes[[i_gene]]$`Pr(>F)`[2]
}
View(mergedDatMarker)
# mergedDatMarker$`BAT(T)_R00_Mean` = rowMeans(mergedDatMarker[,grepl('BAT\\(T\\)\\_R00\\(', colnames(mergedDatMarker))])
# mergedDatMarker$`WAT(T)_R00_Mean` = rowMeans(mergedDatMarker[,grepl('WAT\\(T\\)\\_R00\\(', colnames(mergedDatMarker))])
# mergedDatMarker$`BAT(T)_R20_Mean` = rowMeans(mergedDatMarker[,grepl('BAT\\(T\\)\\_R20\\(', colnames(mergedDatMarker))])
# mergedDatMarker$`iWAT(T)_R20_Mean` = rowMeans(mergedDatMarker[,grepl('iWAT\\(T\\)\\_R20\\(', colnames(mergedDatMarker))])
# mergedDatMarker$BATvsWAT_R00 = mergedDatMarker$`BAT(T)_R00_Mean`/mergedDatMarker$`WAT(T)_R00_Mean`
# mergedDatMarker$BATvsiWAT_R20 = mergedDatMarker$`BAT(T)_R20_Mean`/mergedDatMarker$`iWAT(T)_R20_Mean`

write.xlsx(mergedDatMarker[,c(20:21)], 
           file = paste0(adipocyte_dir,'BrW_Markers/markers_female_male_anova.xls'), 
           sheetName = 'ANOVA', append = FALSE, showNA = FALSE)
write.xlsx(mergedDatMarker, 
           file = paste0(adipocyte_dir,'BrW_Markers/markers_female_male_anova.xls'), 
           sheetName = 'All Information', append = TRUE, showNA = TRUE)


## try one single gene
tt1 = data.frame(ucp1 = dat[,1], tissue = tissue, gender = gender)
summary(tt1)
interaction.plot(tt1$tissue, tt1$gender, tt1$ucp1)
tt = lm(ucp1 ~ tissue + gender, data = tt1)
anova(tt)
summary.aov(aov(ucp1 ~ gender + tissue, data = tt1))


# caculate and display the distance among samples within 1 study
# jiang
# 2017/02/06
## changed 2017/02/06


# ## copy this script from jiang's dir to share/dir
# system('cp /data/home/jiang/projects/scripts/adipocyte/20170126/* /data/home/share/Projects/Adipocyte/rCode/')


library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(birk)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
# code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))

### function
plotFastaQCReadQuality = function(fastaQCFilePath,...){
  dat = read_files_by_line(file_path = fastaQCFilePath)
  skipNumber = grep('>>Per base sequence quality',dat)
  endNumber = grep('>>Per tile sequence quality',dat)-3-skipNumber
  dat= read.delim(file = fastaQCFilePath,
                  sep = '\t',stringsAsFactors = FALSE,row.names = 1,
                  skip = skipNumber,nrows = endNumber)
  datMax = max(dat)
  dat$x = 1:nrow(dat)
  plot(dat$x,dat$Mean,type = 'n',
       ylim = c(0,datMax+0.1),xlim = c(0.5,nrow(dat)+0.5),
       xaxs='i',yaxs='i',col='blue',xaxt='n',
       #main = strsplit(basename(fastaQCFilePath),split = '\\.txt')[[1]][1],
       cex.main=0.8,
       #xlab = 'Position in read (bp)',ylab = 'Phred quality score',
       xlab = '',ylab = '',
       cex.lab=0.5,cex.axis=0.5,...)
  # 
  # axis(side=1,at=seq(1,nrow(dat),2),
  #      labels = rownames(dat)[seq(1,nrow(dat),2)],
  #      las=0,cex.axis=0.5,tck=0)
  axis(side=1,at=1:nrow(dat),line = -1,
       labels = rownames(dat),col = 'white',
       las=0,cex.axis=0.5,tck=0)
  title(ylab="Phred quality score", line=2, cex.lab=0.6)
  title(xlab="Position in read (bp)", line=1, cex.lab=0.8)
  
  ## add background
  for(i in 1:nrow(dat)){
    if(i%%2 == 0){
      rect((dat$x[i]-0.5),28,(dat$x[i]+0.5),datMax+0.1,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Greens")))(255)[230])
      rect((dat$x[i]-0.5),20,(dat$x[i]+0.5),28,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Oranges")))(255)[230])
      rect((dat$x[i]-0.5),0,(dat$x[i]+0.5),20,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Reds")))(255)[230])
      
    }else{
      rect((dat$x[i]-0.5),28,(dat$x[i]+0.5),datMax+0.1,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Greens")))(255)[200])
      rect((dat$x[i]-0.5),20,(dat$x[i]+0.5),28,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Oranges")))(255)[200])
      rect((dat$x[i]-0.5),0,(dat$x[i]+0.5),20,lwd=0,border = NA,
           col =colorRampPalette(rev(brewer.pal(8, "Reds")))(255)[200])
    }
  }
  rect((dat$x-0.4),dat$Lower.Quartile,(dat$x+0.4),dat$Upper.Quartile,
       col = 'yellow',lwd=0.01)
  lines(dat$x,dat$Mean,type = 'l',col='blue')
  segments(dat$x-0.4,dat$X10th.Percentile,dat$x+0.4,dat$X10th.Percentile,col = 'black')
  segments(dat$x-0.4,dat$X90th.Percentile,dat$x+0.4,dat$X90th.Percentile,col = 'black')
  segments(dat$x,dat$Upper.Quartile,dat$x,dat$X90th.Percentile,col = 'black')
  segments(dat$x,dat$Lower.Quartile,dat$x,dat$X10th.Percentile,col = 'black')
  segments(dat$x-0.4,dat$Median,dat$x+0.4,dat$Median,col = 'red')
}


#### plot for Mouse training data, Mouse comfirmation data from Martin and human data used for prediction

fastqcFilesMouse = list.files(paste0(adipocyte_dir,'RNASeq/'),recursive = TRUE,
                         pattern = 'R[0-9]+\\(r[0-9]\\).*.txt',
                         full.names = TRUE)
fastqcFilesMouse = fastqcFilesMouse[!grepl('old|m11R',fastqcFilesMouse)]

##############################
##### training Mouse for tg(after trimming) and fq(Begore trimming)
#############################

fastqcFilesMouseTrain = fastqcFilesMouse[!grepl('m00R|m20R',fastqcFilesMouse)]
#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseTrainingBeforeTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseTrainFastq = fastqcFilesMouseTrain[grepl('\\/fq\\/',fastqcFilesMouseTrain)]
for(i in 1:length(fastqcFilesMouseTrainFastq)){
  plotFastaQCReadQuality(fastqcFilesMouseTrainFastq[i],
                         main = strsplit(basename(fastqcFilesMouseTrainFastq[i]),
                                         split = '\\.txt')[[1]][1])
}
dev.off()


pdf(paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseTrainingAfterTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseTrainTrim = fastqcFilesMouseTrain[grepl('\\/tg\\/',fastqcFilesMouseTrain)]
for(i in 1:length(fastqcFilesMouseTrainTrim)){
  plotFastaQCReadQuality(fastqcFilesMouseTrainTrim[i],
                         main = strsplit(basename(fastqcFilesMouseTrainTrim[i]),
                                         split = '\\.txt')[[1]][1])
}
dev.off()


#########################################################################################################
#################################################
##### for Mouse comfirmation RNAseq, Martin's data both Female and Male mice
###################################################
fastqcFilesMouseConfirmation = fastqcFilesMouse[grepl('m00R|m20R',fastqcFilesMouse)]
#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationBeforeTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationFastq = fastqcFilesMouseConfirmation[grepl('\\/fq\\/',fastqcFilesMouseConfirmation)]
for(i in 1:length(fastqcFilesMouseConfirmationFastq)){
  mainTitle = gsub('_R0','_R00',strsplit(basename(fastqcFilesMouseConfirmationFastq[i]),
                                         split = '\\.txt')[[1]][1])
  mainTitle = gsub('_R20','(T)_R20', mainTitle)
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationFastq[i],
                         main = mainTitle)
}
dev.off()


#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationAfterTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationTrim = fastqcFilesMouseConfirmation[grepl('\\/tg\\/',fastqcFilesMouseConfirmation)]
for(i in 1:length(fastqcFilesMouseConfirmationTrim)){
  mainTitle = gsub('_R0','_R00',strsplit(basename(fastqcFilesMouseConfirmationTrim[i]),
                                         split = '\\.txt')[[1]][1])
  mainTitle = gsub('_R20','(T)_R20', mainTitle)
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationTrim[i],
                         main = mainTitle)
}
dev.off()
#######################################################################################################
#################################################
##### for Mouse comfirmation RNAseq, Martin's data Female mouse
###################################################
fastqcFilesMouseConfirmationFemale = fastqcFilesMouse[grepl('m00R',fastqcFilesMouse)]
#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationBeforeTrimmingFemale.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationFemaleFastq = fastqcFilesMouseConfirmationFemale[grepl('\\/fq\\/',fastqcFilesMouseConfirmationFemale)]
for(i in 1:length(fastqcFilesMouseConfirmationFemaleFastq)){
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationFemaleFastq[i],
                         main = gsub('_R0','_R00',strsplit(basename(fastqcFilesMouseConfirmationFemaleFastq[i]),
                                         split = '\\.txt')[[1]][1]))
}
dev.off()


#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationAfterTrimmingFemale.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationFemaleTrim = fastqcFilesMouseConfirmationFemale[grepl('\\/tg\\/',fastqcFilesMouseConfirmationFemale)]
for(i in 1:length(fastqcFilesMouseConfirmationFemaleTrim)){
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationFemaleTrim[i],
                         main = gsub('_R0','_R00',strsplit(basename(fastqcFilesMouseConfirmationFemaleTrim[i]),
                                         split = '\\.txt')[[1]][1]))
}
dev.off()
#################################################
#################################################
##### for Male Mouse comfirmation RNAseq, Martin's data 
###################################################
fastqcFilesMouseConfirmationMale = fastqcFilesMouse[grepl('m20R',fastqcFilesMouse)]
#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationBeforeTrimmingMale.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationMaleFastq = fastqcFilesMouseConfirmationMale[grepl('\\/fq\\/',fastqcFilesMouseConfirmationMale)]
for(i in 1:length(fastqcFilesMouseConfirmationMaleFastq)){
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationMaleFastq[i],
                         main = gsub('_R','(T)_R',strsplit(basename(fastqcFilesMouseConfirmationMaleFastq[i]),
                                                         split = '\\.txt')[[1]][1]))
}
dev.off()



#pdf(paste0(plots_dir,'RNAseqReadsQuality.pdf'),
pdf(paste0(adipocyte_dir,
           'RNASeq/All_BeBrW/tg/RNAseqReadsQualityMouseConfirmationAfterTrimmingMale.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesMouseConfirmationMaleTrim = fastqcFilesMouseConfirmationMale[grepl('\\/tg\\/',fastqcFilesMouseConfirmationMale)]
for(i in 1:length(fastqcFilesMouseConfirmationMaleTrim)){
  plotFastaQCReadQuality(fastqcFilesMouseConfirmationMaleTrim[i],
                         main = gsub('_R','(T)_R',strsplit(basename(fastqcFilesMouseConfirmationMaleTrim[i]),
                                                         split = '\\.txt')[[1]][1]))
}
dev.off()

###############################################
### human data
##############################################3
fastqcFilesHuman = list.files(paste0(adipocyte_dir,'RNASeq/'),recursive = TRUE,
                              pattern = 'hR.*.txt',
                              full.names = TRUE)
fastqcFilesHuman = fastqcFilesHuman[!grepl('old',fastqcFilesHuman)]

pdf(paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/RNAseqReadsQualityHumanBeforeTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesHumanFastq = fastqcFilesHuman[grepl('\\/fq\\/',fastqcFilesHuman)]
for(i in 1:length(fastqcFilesHumanFastq)){
  plotFastaQCReadQuality(fastqcFilesHumanFastq[i],
                         main = strsplit(basename(fastqcFilesHumanFastq[i]),
                                         split = '\\.txt')[[1]][1])
}
dev.off()

pdf(paste0(adipocyte_dir,'RNASeq/All_BeBrW/tg/RNAseqReadsQualityHumanAfterTrimming.pdf'),
    family="Helvetica",width = 18.3/2.54,height = 24.7/2.54)
#par(mfrow=c(3,2))
par(mfrow=c(4,2),mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
fastqcFilesHumanTrim = fastqcFilesHuman[grepl('\\/tg\\/',fastqcFilesHuman)]
for(i in 1:length(fastqcFilesHumanTrim)){
  plotFastaQCReadQuality(fastqcFilesHumanTrim[i],
                         main = strsplit(basename(fastqcFilesHumanTrim[i]),
                                         split = '\\.txt')[[1]][1])
}
dev.off()

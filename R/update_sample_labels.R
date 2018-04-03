# update sample labels
# jiang
# 2017/09/28


library(xlsx)
library(hashmap)
library(stringi)
library(stringr)

 #adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
# code_path = '/data/home/share/Projects/Adipocyte/rCode/'

adipocyte_dir = '//home/jiang/MitoShare/Adipocyte/'


#########################################################

##########################################################
### get old labels
# #human 
# human_labels_old = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/TableS2_Cheng_old_label.xls'), 
#                     sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
# human_labels_old = human_labels_old[!is.na(human_labels_old$Study.ID), c("Sample.ID", "Sample.ID.1")]
# colnames(human_labels_old) = c('sample','label_old')
# 
# # mouse
# # mouse_labels_old = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/TableS1_Cheng_old_label.xls'), 
# #                              sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
# mouse_labels_old = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/TableS1_old.xls'), 
#                              sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
# mouse_labels_old = mouse_labels_old[!is.na(mouse_labels_old$Study.ID), c( "Sample.ID", "Study.ID")]
# colnames(mouse_labels_old) = c('sample','label_old')

######################
label_old_part = read.delim(file =  paste0(adipocyte_dir, 'DataInfo/sampleLabels/All_Dataset.txt'), 
                            header = TRUE, stringsAsFactors = FALSE)
label_old_part = label_old_part[, c( "SampleID", "StudyID")]
colnames(label_old_part) = c('sample','label_old')
# mouse M13
m13_labels_old = read.table(file =  paste0(adipocyte_dir, 'DataInfo/MicroArray_Dataset.NoCore.txt'), 
                            header = TRUE, stringsAsFactors = FALSE)
m13_sample = strsplit(m13_labels_old$Samples, split = ';')[[1]]
m13_labels = strsplit(m13_labels_old$TissueName, split = ';')[[1]]
m13_labels_old = data.frame(sample = m13_sample, label_old = m13_labels)

#labels_old = rbind(mouse_labels_old, m13_labels_old, human_labels_old)
#rm(mouse_labels_old, m13_labels_old, human_labels_old)

labels_old = rbind(label_old_part, m13_labels_old)
#############################################################
#### get new labels
#human 
human_labels_new = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/TableS2.xls'), 
                             sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
human_labels_new = human_labels_new[!is.na(human_labels_new$Study.ID), c("Sample.ID", "Dataset.ID")]
colnames(human_labels_new) = c('sample','label_new')

# mouse
mouse_labels_new = read.xlsx(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/TableS1_20171002.xls'), 
                             sheetIndex = 1, stringsAsFactors = FALSE, header = TRUE)
mouse_labels_new = mouse_labels_new[!is.na(mouse_labels_new$Study.ID), c( "Sample.ID", "Dataset.ID")]
colnames(mouse_labels_new) = c('sample','label_new')

labels_new = rbind(mouse_labels_new, human_labels_new)
rm(mouse_labels_new, human_labels_new)

labels_map = merge(labels_old, labels_new, by = 'sample')
dim(labels_map)
length(unique(labels_map$sample))

### change M1 to M01
change_format = function(tmp_string){
    if(grepl('hM[0-9]$|_M[0-9][^0-9]',tmp_string)){
        tmp = str_split(tmp_string, pattern = 'M')[[1]]
        return(paste0(tmp[1],'M0',tmp[2]))
    }else{
        return(tmp_string)
    }
}

labels_map$label_old2 = sapply(labels_map$label_old, function(x)change_format(x))
write.table(labels_map, file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/samples_labels_20171002.txt'),
            sep = '\t', quote = TRUE, row.names = FALSE)

###################################################
# the 1st time update

sample_label_map = hashmap(labels_map$label_old2, labels_map$label_new)
#####################################################################
# # the second time update
# 
oldLabels = read.table(file = paste0(adipocyte_dir, 'DataInfo/sampleLabels/samples_labels.txt'),
                       sep = '\t', header = TRUE)
labels_map_20171002 = merge(oldLabels, labels_map[,c(1,3)], by = 'sample')
sample_label_map = hashmap(labels_map_20171002$label_new.x, labels_map_20171002$label_new.y)


#####################################################################################
## update sample labels
update_labels = function(filePath, dictionary){
  tmp_file = read.table(file = filePath, header = TRUE, check.names = FALSE, sep = '\t')
  if(grepl('pca$', basename(filePath))){
      if(sum(is.na(dictionary[[rownames(tmp_file)]]))>0){
          print(basename(filePath))
          print(rownames(tmp_file)[1:3])
      }
  }else{
      if(sum(is.na(dictionary[[colnames(tmp_file)]]))>0){
          print(basename(filePath))
          print(colnames(tmp_file)[1:3])
      }else{
          colnames(tmp_file) = dictionary[[colnames(tmp_file)]]
          write.table(tmp_file,file = filePath,
                      sep = '\t')
      }
  }
}


file_folder = '/home/jiang/Downloads/adipose/data/'

files = list.files(path = paste0(file_folder, 'MicroArray/'), 
                   full.names = TRUE, pattern = '.txt$')
length(files)
update_labels(filePath = files[1], dictionary = sample_label_map)
for(file in files){
    update_labels(filePath = file, dictionary = sample_label_map)
}

files = list.files(path = paste0(file_folder, 'RNASeq/'), 
                   full.names = TRUE, pattern = '.txt$')
length(files)

for(file in files){
    update_labels(filePath = file, dictionary = sample_label_map)
}




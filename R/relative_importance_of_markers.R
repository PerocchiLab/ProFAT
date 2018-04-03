library(stringr)
library(caret)
library(hashmap)
library(sva)

project_dir<-'/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
source(paste0(code_path,'gar_fun.r'))
source(paste0(project_dir, "rCode/Functions_Cheng.R"))


###############################################

#######################################
## neural network model
#######################################
# Different machine learning algorithms.
#ml_method<-c("rf", "nnet", "nb", "glm", "rpart", "svmRadialSigma")
pred_method = 'nnet'
pre_process<-NULL # NULL, scale, center. corr
t_length<-3 # grid tune length.

########## read data
raw_brw_file<-paste0(project_dir, "/BrW_Markers/data_BrW.txt")
raw_train_data<-read.table(raw_brw_file, sep="\t", check.names=FALSE, row.names=1)
batch_VM = unname(sapply(colnames(raw_train_data), function(x) str_extract_all(x,'M[0-9]+|R[0-9]+')[[1]]))
train_data_normalized = ComBat(dat = raw_train_data, batch = batch_VM)

# Load the markers ...
marker_dir<-paste0(project_dir, "BrW_Markers")
marker_file<-paste0(marker_dir, "/marker.txt")
marker_data<-read.table(marker_file, sep="\t")
marker_data<-marker_data[order(marker_data$logFC), ]

raw_train_data_marker = raw_train_data[rownames(raw_train_data) %in% rownames(marker_data),]
train_data_normalized_marker = train_data_normalized[rownames(train_data_normalized) %in% rownames(marker_data),]


#train_data = data.frame(t(raw_train_data_marker))
train_data = data.frame(t(train_data_normalized_marker))
# train_data = raw_train_data_marker[,c(1:50)]
# train_data = data.frame(t(train_data))
# test_dat = raw_train_data_marker[,c(51:59)]
# test_dat = data.frame(t(test_dat))

tissue_type<-str_extract(rownames(train_data), "BAT|WAT")
#train_data<-as.data.frame(t(train_data))
train_data$TissueType<-tissue_type

########################################################
## train the model
control <- trainControl(method="LOOCV", classProbs=TRUE)
model <- train(TissueType~., data=train_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength = t_length)

#############################################################
# plot variable importance

# get variable importance
weights = data.frame(gar.fun(train_data$TissueType,model,bar.plot = FALSE))
relative_weight = weights[, c("inp.cont", "rel.imp")]
relative_weight = relative_weight[order(relative_weight$rel.imp, decreasing = FALSE),]
relative_weight$features = factor(row.names(relative_weight)[order(relative_weight$rel.imp, decreasing = FALSE)], levels = row.names(relative_weight)[order(relative_weight$rel.imp, decreasing = FALSE)])


## set colors according to ENSEMBLE ID
markerFolder = paste0(project_dir,'BrW_Markers/')
## get colors for markers
mouseMarker = read.delim(file = paste0(markerFolder,'marker.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Ensembl = rownames(mouseMarker)
mouseMarkerBat = read.delim(file = paste0(markerFolder,'marker.BAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarkerWat = read.delim(file = paste0(markerFolder,'marker.WAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Type = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'BAT','WAT'))
mouseMarker$Color = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'brown2','papayawhip'))
hashmap_type = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$Type)
hashmap_colors = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$Color)
hashmap_names = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$GeneName)

relative_weight$type = factor(hashmap_type[[relative_weight$features]])
relative_weight$colors = hashmap_colors[[relative_weight$features]]
relative_weight$gene = hashmap_names[[relative_weight$features]]
relative_weight$gene = factor(relative_weight$gene[order(relative_weight$rel.imp, decreasing = FALSE)], levels = relative_weight$gene[order(relative_weight$rel.imp, decreasing = FALSE)])


pdf(paste0(markerFolder, 'relative_importance_markers_using_normalized_data.pdf'))

ggplot(relative_weight, aes(x = gene, y = rel.imp, fill = type)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('brown','papayawhip')) +
  theme(legend.position="none") + # remove legend 
  theme(axis.text.x = element_text(size = 10, angle = 45)) +
  coord_flip(ylim = c(-1,1)) +
  ylab('Relative Importance') +
  xlab('Marker Genes')
dev.off()

# Predict the test data
pred<-predict(model, test_dat, type="raw")


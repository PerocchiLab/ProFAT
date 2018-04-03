# Predict the datasets
library(caret)
library(stringr)
library(sva)
library(doMC)
library(biomaRt)
library(pheatmap)
library(e1071)
library(nnet)
library(randomForest)
library(rpart)
library(assertthat)
library(neuralnet)
library(gplots)
library(amap)

project_dir<-"/data/home/cheng/Adipocyte/"

# Input parameters
platform<-"MicroArray"  # data platform, one of "MicroArray", "RNASeq"
ml_method<-"nnet" # Machine learning algorithm. One of "rf", "nnet", "nb", "glm", "rpart", "svmRadialSigma"

########
# test study name, as an example: h11M, it is composed of three parts:
# h-human,m-mouse; 11-study number; R-RNAseq,M-MicroArray
test_study<-"h11M"
########

################################################################################
# Predict the data together with Training Data.
################################################################################
# Load the Raw training data ...
raw_brw_file<-paste0(project_dir, "BrW_Markers/data_BrW.txt")
raw_train_data<-read.table(raw_brw_file, sep="\t", check.names=FALSE, row.names=1)

# Load the markers ...
marker_file<-paste0(project_dir, "BrW_Markers/marker.txt")
marker_data<-read.table(marker_file, sep="\t")
marker_data<-marker_data[order(marker_data$logFC), ]


study_spe<-"mouse"
if (grepl("^h", test_study)){
    study_spe<-"human"
}
study_dir<-paste0(project_dir, platform, "/", test_study)

# Load each individual study
if (platform == "MicroArray"){
    test_data_file<-paste0(study_dir, "/ND_", test_study, ".txt")
} else if (platform == "RNASeq") {
    test_data_file<-paste0(study_dir, "/rlog.txt")
}
test_pre_out<-paste0(study_dir, "/pred.", ml_method, ".txt")
test_data<-read.table(test_data_file, sep="\t", check.names=FALSE, row.names=1)

if (study_spe == "human"){
    ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
    gene_name = getBM(attributes=c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"), filters="ensembl_gene_id",
                      values=row.names(test_data), mart=ensembl_mart)
    gene_name<-subset(gene_name, mmusculus_homolog_orthology_type=="ortholog_one2one")
    rownames(gene_name)<-gene_name$ensembl_gene_id
    test_data<-test_data[gene_name$ensembl_gene_id, ]
    row.names(test_data)<-gene_name[rownames(test_data), "mmusculus_homolog_ensembl_gene"]
}
cnames_test_data<-colnames(test_data)

# Merge with Raw training data.
TT_data<-merge(test_data, raw_train_data, by=0)
rownames(TT_data)<-TT_data$Row.names
TT_data<-subset(TT_data, select=-Row.names)

batch_VM<-gsub(".*_([MR]\\d+)\\(r\\d+\\)", "\\1", colnames(TT_data))
batch_VM<-gsub(".*_", "", batch_VM)

if (length(unique(batch_VM))>1){
    batch_TT_data<-ComBat(dat = TT_data, batch = batch_VM)
} else {
    batch_TT_data<-TT_data
}

# Use the markers that are both in Test and Train data
# markers_used<-intersect(rownames(train_data_all), rownames(batch_test_data))
markers_used<-intersect(rownames(marker_data), rownames(batch_TT_data))

# Prepare the testing data
batch_test_data<-batch_TT_data[markers_used, cnames_test_data]
batch_test_data<-as.data.frame(t(batch_test_data))

# Prepare the training data
batch_train_data<-batch_TT_data[markers_used, colnames(raw_train_data)]
batch_train_data<-as.data.frame(t(batch_train_data))

tissue_type<-str_extract(rownames(batch_train_data), "BAT|WAT")
batch_train_data$TissueType<-tissue_type

# Train the Model
control <- trainControl(method="LOOCV", classProbs=TRUE)
model <- train(TissueType~., data=batch_train_data, method=ml_method, preProcess=NULL, trControl=control, tuneLength = 3)

# Predict the test data
pred<-predict(model, batch_test_data, type="prob")
write.table(pred, test_pre_out, sep="\t", quote=FALSE)

# Predict the datasets
library(caret)
library(stringr)
library(sva)
library(doMC)
library(biomaRt)
library("pheatmap")
library(e1071)
library(nnet)
library(randomForest)
library(rpart)
library(assertthat)
library(neuralnet)
library(gplots)
library(amap)
library(pROC)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))
fig_dir<-paste0(project_dir, "MLResults/")

# Different data technology
platforms<-c("RNASeq", "MicroArray")

# Different machine learning algorithms.
# gbm: The dataset size is too small
ml_method<-c("rf", "nnet", "nb", "glm", "rpart", "svmRadialSigma")

# data preprocessing ...
pre_process<-NULL # NULL, scale, center. corr
if (is.null(pre_process)){
    pre_pro_txt<-"null"
} else {
    pre_pro_txt<-pre_process
}

type_analysis<-"train" # train
if (type_analysis=="all"){
    studies_m<-c("m01M_BrW", "m02M_BrW", "m03M_BeW", "m04M_BrW", "m05M_BrW", "m06M_BrW", "m07M_BrW", "m08M_BeW", "m09M_BrW", "m10M_BeBrW", "m11M_BeBr", "m12M_BrW", "m13M_BeBrW") #
    studies_r<-c("m01R_BeBr", "m02R_BeBrW", "m03R_BeBrW", "m04R_BrW", "m00R_BrBeW",  "m20R_BrW")
} else if (type_analysis=="train") {
    studies_m<-c("m03M_BeW", "m08M_BeW", "m09M_BrW", "m10M_BeBrW", "m11M_BeBr", "m13M_BeBrW")
    studies_r<-c("m01R_BeBr", "m02R_BeBrW", "m03R_BeBrW", "m00R_BrBeW",  "m20R_BrW")
}

sensi<-data.frame(ML_Method=as.character(NULL), Accu=as.numeric(NULL), Dataset=as.character(NULL))

treated_as_BAT<-0
data_auc_nnet<-data.frame(Response=as.character(NULL), Prediction=as.numeric(NULL))

for (pred_method in ml_method){
    for (t_study in c(studies_r, studies_m)){
        # cat(pred_method, "\t", t_study, "\n")
        if (grepl("R_", t_study)){
            data_platform<-"RNASeq"
            pre_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/all.pred.wt.null.", pred_method, ".txt")
        } else if (grepl("M_", t_study)){
            data_platform<-"MicroArray"
            pre_file<-paste0(project_dir, data_platform, "/", t_study, "/all.pred.wt.null.", pred_method, ".txt")
        }
        pre_data<-read.table(pre_file, sep="\t", header=TRUE)

        if (treated_as_BAT == 0){
            pre_data<-pre_data[!grepl("WAT\\(", rownames(pre_data)), ]
        }

        # Only retained the data for prediction...
        pre_data<-pre_data[complete.cases(pre_data), ]
        study_str<-get_study_str(t_study)
        pre_data<-pre_data[grepl(paste0(study_str, "\\("), rownames(pre_data)), ]

        total_samples<-dim(pre_data)[1]

        # the true tissue type
        pre_data$Tissue_Type<-"NA"
        pre_data[grepl("WAT", rownames(pre_data)), "Tissue_Type"]<-"WAT"
        pre_data[grepl("BAT", rownames(pre_data)), "Tissue_Type"]<-"BAT"
        # treated as BAT
        if (treated_as_BAT == 1){
            pre_data[grepl("WAT\\(", rownames(pre_data)), "Tissue_Type"]<-"BAT"
        }

        # the predicted tissue type
        pre_data$Pre_Type<-"NA"
        pre_data[pre_data$BAT>0.5, "Pre_Type"]<-"BAT"
        pre_data[pre_data$BAT<0.5, "Pre_Type"]<-"WAT"

        # filter out the NAs
        pre_data<-subset(pre_data, Tissue_Type != "NA" & Pre_Type !="NA")
        study_accu<- sum(pre_data$Tissue_Type==pre_data$Pre_Type)/total_samples
        study_accu_df<-data.frame(ML_Method=pred_method, Accu=study_accu, Dataset=t_study)

        sensi<-rbind(sensi, study_accu_df)

        # data for AUC
        if (pred_method == "nnet"){
            pre_data_auc<-pre_data[, c("Tissue_Type", "BAT")]
            colnames(pre_data_auc)<-c("Response", "Prediction")
            data_auc_nnet<-rbind(data_auc_nnet, pre_data_auc)
        }
    }
}

if (treated_as_BAT == 1){
    fig_out<-paste0(fig_dir, "algorithms.wt.beige.", type_analysis, ".eps")
    res_out<-paste0(fig_dir, "algorithms.wt.beige.", type_analysis, ".txt")

} else {
    fig_out<-paste0(fig_dir, "algorithms.rm.beige.", type_analysis, ".eps")
    res_out<-paste0(fig_dir, "algorithms.rm.beige.", type_analysis, ".txt")
}

setEPS()
postscript(fig_out)

levels(sensi$ML_Method)<-gsub("RadialSigma", "", levels(sensi$ML_Method))
boxplot(Accu ~ ML_Method, data=sensi, ylim=c(0,1))
dev.off()

write.table(sensi, res_out, sep="\t", quote=FALSE, row.names=FALSE)

# Generate AUC
roc_obj <- roc(data_auc_nnet$Response, data_auc_nnet$Prediction)

if (treated_as_BAT == 1){
    fig_out<-paste0(fig_dir, "AUC.nnet.wt.beige.", type_analysis, ".eps")
} else {
    fig_out<-paste0(fig_dir, "AUC.nnet.rm.beige.", type_analysis, ".eps")
}
setEPS()
postscript(fig_out)
plot(roc_obj$sensitivities, roc_obj$specificities, type="l", xlab="Sensitivity", ylab="Specificity")
text(0.5, 0.5, paste0("AUC is ", round(roc_obj$auc, 4)))
dev.off()


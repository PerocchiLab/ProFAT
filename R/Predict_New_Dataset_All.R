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
library(gPCA)
library(amap)

registerDoMC(cores = 16)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))
platforms<-c("RNASeq", "MicroArray")

# Different machine learning algorithms.
ml_method<-c("rf", "nnet", "nb", "glm", "rpart", "svmRadialSigma")
# ml_method<-c("nnet")
#ml_method<-c("svmRadialSigma")

# gbm: The dataset size is too small
pre_process<-NULL # NULL, scale, center. corr
t_length<-3 # grid tune length.

if (is.null(pre_process)){
    pre_pro_txt<-"null"
} else {
    pre_pro_txt<-pre_process
}

################################################################################
# Predict the data together with Training Data.
################################################################################
# Load the training data ...
brw_file<-paste0(project_dir, "/BrW_Markers/rbatch.data_BrW.txt")
train_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)
raw_brw_file<-paste0(project_dir, "/BrW_Markers/data_BrW.txt")
raw_train_data<-read.table(raw_brw_file, sep="\t", check.names=FALSE, row.names=1)

# Load the markers ...
marker_dir<-paste0(project_dir, "BrW_Markers")
marker_file<-paste0(marker_dir, "/marker.txt")
marker_data<-read.table(marker_file, sep="\t")
marker_data<-marker_data[order(marker_data$logFC), ]
for (pred_method in ml_method){
    for (data_platform in platforms){
        # Train and Test each study
        test_studies<-dir(paste0(project_dir, data_platform), pattern="^[mh]\\d+[MR]*", full.names=FALSE)
        for (t_study in test_studies){
            if (t_study != "m13M_BeBrW" && t_study!="m20R_BrW"){
                next
            }
            cat(paste("**WithTrainingData******", pred_method, data_platform, t_study, "\n"))
            study_spe<-"mouse"
            if (grepl("^h", t_study)){
                study_spe<-"human"
            }
            # Load each individual study
            if (data_platform == "MicroArray"){
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/ND_", t_study, ".txt")
                test_data_with_train_file<-paste0(project_dir, data_platform, "/", t_study, "/ND_", t_study, ".with.training.txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/all.pred.wt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/all.pheatmap.wt.pdf")
                test_agg_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/all.pheatmap.wt.agg.pdf")
                test_agg_1_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/all.pheatmap.wt.agg1.pdf")
                test_agg_2_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/heatmap.eps")
            } else {
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/rlog.txt")
                test_data_with_train_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/rlog.with.training.txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/all.pred.wt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/all.pheatmap.wt.pdf")
                test_agg_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/all.pheatmap.wt.agg.pdf")
                test_agg_1_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/all.pheatmap.wt.agg1.pdf")
                test_agg_2_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/heatmap.eps")
            }

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
            col_names_study<-colnames(test_data)

            # Merge with Raw training data.
            if (all(colnames(test_data) %in% colnames(raw_train_data))){
                test_data<-raw_train_data
            } else {
                test_data<-merge(test_data, raw_train_data, by=0)
                rownames(test_data)<-test_data$Row.names
                test_data<-subset(test_data, select=-Row.names)
            }

            batch_VM<-gsub(".*_([MR]\\d+)\\(r\\d+\\)", "\\1", colnames(test_data))
            batch_VM<-gsub(".*_", "", batch_VM)

            if (length(unique(batch_VM))>1){
                batch_test_data<-ComBat(dat = test_data, batch = batch_VM)
            } else {
                batch_test_data<-test_data
            }

            write.table(batch_test_data, test_data_with_train_file, sep="\t", quote=FALSE)

            # Use the markers that are both in Test and Train data
            # markers_used<-intersect(rownames(train_data_all), rownames(batch_test_data))
            markers_used<-intersect(rownames(marker_data), rownames(batch_test_data))

            # Prepare the testing data
            batch_test_data<-batch_test_data[markers_used, ]
            batch_test_data<-as.data.frame(t(batch_test_data))

            # Prepare the training data
            #train_data<-train_data_all[rownames(marker_data), ]
            train_data<-batch_test_data[!rownames(batch_test_data) %in% col_names_study, ]

            tissue_type<-str_extract(rownames(train_data), "BAT|WAT")
            #train_data<-as.data.frame(t(train_data))
            train_data$TissueType<-tissue_type

            # Train the Model
            if (!file.exists(test_pre_out)){
                control <- trainControl(method="LOOCV", classProbs=TRUE)
                model <- train(TissueType~., data=train_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength = t_length)

                # Predict the test data
                pred<-predict(model, batch_test_data, type="prob")
                if (!identical(rownames(pred), rownames(batch_test_data))){
                    rownames(pred)<-rownames(batch_test_data)
                }
                write.table(pred, test_pre_out, sep="\t", quote=FALSE)
            }

            #########################
            # Plot the testing data.
            #########################
            # Plot all the samples seperately
            #pheatmap_plot(batch_test_data, test_fig_out)

            # Plot all the aggregated data together.
            batch_test_data$Samples<-gsub("\\(r\\d+\\)$", "", rownames(batch_test_data))
            batch_test_data_agg<-aggregate(.~Samples, data=batch_test_data, mean)
            row.names(batch_test_data_agg)<-batch_test_data_agg[,1]
            batch_test_data_agg<-batch_test_data_agg[,-1]
            #pheatmap_plot(batch_test_data_agg, test_agg_fig_out)

            #Aggreate all the training data into two values and keep all the samples
            batch_test_data<-subset(batch_test_data, select=-Samples)
            s_data<-batch_test_data[rownames(batch_test_data) %in% col_names_study, ]
            train_data_agg<-aggregate(.~TissueType, data=train_data, mean)
            row.names(train_data_agg)<-train_data_agg[,1]
            train_data_agg<-train_data_agg[,-1]
            #pheatmap_plot(rbind(train_data_agg, s_data), test_agg_1_fig_out, use_heatmap=0)

            # Plot the data using Heatmap function...
            pheatmap_plot(rbind(train_data_agg, s_data), test_agg_2_fig_out, use_heatmap=1)
        }
    }
}


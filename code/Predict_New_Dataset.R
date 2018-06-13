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

registerDoMC(cores = 16)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

marker_dir<-paste0(project_dir, "MarkerPrediction")

platforms<-c("RNASeq", "MicroArray")
#platforms<-c("RNASeq")

marker_file<-paste0(marker_dir, "/marker.txt")
marker_data<-read.table(marker_file, sep="\t")
t_length<-3 # grid tune length.

# Different machine learning algorithms.
ml_method<-c("rf", "nnet", "nb", "glm", "rpart", "svmRadialSigma")
#ml_method<-c("svmRadialSigma")

# gbm: The dataset size is too small
pre_process<-"scale" # NULL, scale, center. corr

if (is.null(pre_process)){
    pre_pro_txt<-"null"
} else {
    pre_pro_txt<-pre_process
}

################################################################################
# Predict the data together with Training Data.
################################################################################
for (pred_method in ml_method){
    # Test each individual study
    for (data_platform in platforms){
        brw_file<-paste0(project_dir, data_platform, "/BrW_Markers/rbatch.data_BrW.txt")
        train_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)

        raw_brw_file<-paste0(project_dir, data_platform, "/BrW_Markers/data_BrW.txt")
        raw_train_data<-read.table(raw_brw_file, sep="\t", check.names=FALSE, row.names=1)

        # Train and Test each study
        test_studies<-dir(paste0(project_dir, data_platform), pattern="^[mh]\\d+[MR]*", full.names=FALSE)
        for (t_study in test_studies){
            cat(paste("**WithTrainingData******", pred_method, data_platform, t_study, "\n"))
            study_spe<-"mouse"
            if (grepl("^h", t_study)){
                study_spe<-"human"
            }
            # Load each individual study
            if (data_platform == "MicroArray"){
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/ND_", t_study, ".txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/pred.wt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/pheatmap.wt.pdf")
            } else {
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/rlog.txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/pred.wt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/pheatmap.wt.pdf")
            }

            if (file.exists(test_pre_out)){
                next
            }

            test_data<-read.table(test_data_file, sep="\t", check.names=FALSE, row.names=1)
            if (study_spe == "human"){
                ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
                gene_name = getBM(attributes=c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"), filters="ensembl_gene_id",
                                  values=row.names(test_data), mart=ensembl_mart)
                gene_name<-subset(gene_name, mmusculus_homolog_orthology_type=="ortholog_one2one")
                test_data<-test_data[gene_name$ensembl_gene_id, ]
                row.names(test_data)<-gene_name$mmusculus_homolog_ensembl_gene
            }

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

            # Use the markers that are both in Test and Train data
            markers_used<-intersect(rownames(train_data_all), rownames(batch_test_data))
            markers_used<-intersect(rownames(marker_data), markers_used)

            # Prepare the testing data
            batch_test_data<-batch_test_data[markers_used, ]
            batch_test_data<-as.data.frame(t(batch_test_data))

            # Prepare the training data
            train_data<-train_data_all[markers_used, ]
            tissue_type<-str_extract(colnames(train_data), "BAT|WAT")
            train_data<-as.data.frame(t(train_data))
            train_data$TissueType<-tissue_type

            # Train the Model
            control <- trainControl(method="LOOCV", classProbs=TRUE)
            model <- train(TissueType~., data=train_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength = t_length)

            # Predict the test data
            pred<-predict(model, batch_test_data, type="prob")
            if (!identical(rownames(pred), rownames(batch_test_data))){
                rownames(pred)<-rownames(batch_test_data)
            }
            write.table(pred, test_pre_out, sep="\t", quote=FALSE)

            #########################
            # Plot the testing data.
            #########################
            pheatmap_plot(batch_test_data, test_fig_out)
        }
    }
}

################################################################################
# Predict the data All_BeBrW
################################################################################
for (pred_method in ml_method){
    # Test All_BeBrW together
    # Test each individual study
    for (data_platform in platforms){
        cat(paste("**All_BeBrW******", pred_method, data_platform, "\n"))

        brw_file<-paste0(project_dir, data_platform, "/BrW_Markers/rbatch.data_BrW.txt")
        train_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)

        if (data_platform == "MicroArray"){
            test_data_file<-paste0(project_dir, data_platform, "/All_BeBrW/rbatch.All.ND.txt")
            test_pre_out<-paste0(project_dir, data_platform, "/All_BeBrW/pred.", pre_pro_txt, ".", pred_method, ".txt")
        } else {
            test_data_file<-paste0(project_dir, data_platform, "/All_BeBrW/tg/rlog.txt")
            test_pre_out<-paste0(project_dir, data_platform, "/All_BeBrW/tg/pred.", pre_pro_txt, ".", pred_method, ".txt")
        }

        if (file.exists(test_pre_out)){
            next
        }

        test_data<-read.table(test_data_file, sep="\t", check.names=FALSE, row.names=1)

        # Use the markers that are both in Test and Train data
        markers_used<-intersect(rownames(train_data_all), rownames(test_data))
        markers_used<-intersect(rownames(marker_data), markers_used)

        # Prepare the testing data
        test_data<-test_data[markers_used, ]
        test_data<-as.data.frame(t(test_data))

        # Prepare the training data
        train_data<-train_data_all[markers_used, ]
        tissue_type<-str_extract(colnames(train_data), "BAT|WAT")
        train_data<-as.data.frame(t(train_data))
        train_data$TissueType<-tissue_type

        # Train the Model
        control <- trainControl(method="LOOCV", classProbs=TRUE)
        model <- train(TissueType~., data=train_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength=t_length)

        # Predict the test data
        pred<-predict(model, test_data, type="prob")
        if (!identical(rownames(pred), rownames(test_data))){
            rownames(pred)<-rownames(test_data)
        }

        write.table(pred, test_pre_out, sep="\t", quote=FALSE)
    }
}

################################################################################
# Predict the test data ALONE, NOT together with training data.
################################################################################
for (pred_method in ml_method){
    # Test each individual study
    for (data_platform in platforms){
        brw_file<-paste0(project_dir, data_platform, "/BrW_Markers/rbatch.data_BrW.txt")
        train_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)

        # Train and Test each study
        test_studies<-dir(paste0(project_dir, data_platform), pattern="^[mh]\\d+[MR]*", full.names=FALSE)
        for (t_study in test_studies){
            cat(paste("**WithoutTrainingData******", pred_method, data_platform, t_study, "\n"))
            study_spe<-"mouse"
            if (grepl("^h", t_study)){
                study_spe<-"human"
            }
            if (data_platform == "MicroArray"){
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/ND_", t_study, ".txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/pred.nt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/pheatmap.pdf")
            } else {
                test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/rlog.txt")
                test_pre_out<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/pred.nt.", pre_pro_txt, ".", pred_method, ".txt")
                test_fig_out<-paste0(project_dir, data_platform, "/", t_study, "/tg/pheatmap.pdf")
            }

            if (file.exists(test_pre_out)){
                next
            }

            test_data<-read.table(test_data_file, sep="\t", check.names=FALSE, row.names=1)

            if (study_spe == "human"){
                ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
                gene_name = getBM(attributes=c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"), filters="ensembl_gene_id",
                                  values=row.names(test_data), mart=ensembl_mart)
                gene_name<-subset(gene_name, mmusculus_homolog_orthology_type=="ortholog_one2one")
                test_data<-test_data[gene_name$ensembl_gene_id, ]
                row.names(test_data)<-gene_name$mmusculus_homolog_ensembl_gene
            }

            # Use the markers that are both in Test and Train data
            markers_used<-intersect(rownames(train_data_all), rownames(test_data))
            markers_used<-intersect(rownames(marker_data), markers_used)

            # Prepare the testing data
            test_data<-test_data[markers_used, ]
            test_data<-as.data.frame(t(test_data))

            # Prepare the training data
            train_data<-train_data_all[markers_used, ]
            tissue_type<-str_extract(colnames(train_data), "BAT|WAT")
            train_data<-as.data.frame(t(train_data))

            # Plot the heatmap of the training and testing data.
            pheatmap_plot(test_data, test_fig_out, train_data, paste0(project_dir, data_platform, "/BrW_Markers/pheatmap.pdf"))

            train_data$TissueType<-tissue_type
            # Train the Model
            control <- trainControl(method="LOOCV", classProbs=TRUE)
            model <- train(TissueType~., data=train_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength=t_length)

            # Predict the test data
            pred<-predict(model, test_data, type="prob")
            if (!identical(rownames(pred), rownames(batch_test_data))){
                rownames(pred)<-rownames(test_data)
            }

            write.table(pred, test_pre_out, sep="\t", quote=FALSE)
        }
    }
}

################################################################################
# Predict the test data Using the machine learning algorithms Seperately.
################################################################################
# Test each individual study
# for (data_platform in platforms){
#     brw_file<-paste0(project_dir, data_platform, "/BrW_Markers/rbatch.data_BrW.txt")
#     train_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)
#
#     # Train and Test each study
#     test_studies<-dir(paste0(project_dir, data_platform), pattern="^[mh]\\d+[MR]*", full.names=FALSE)
#     for (t_study in test_studies){
#         if (t_study == "m00R_BrBeW"){
#             next
#         }
#         cat(paste("**WithoutTrainingData*****MLSeperate*****", data_platform, t_study, "\n"))
#         study_spe<-"mouse"
#         if (grepl("^h", t_study)){
#             study_spe<-"human"
#         }
#         if (data_platform == "MicroArray"){
#             test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/ND_", t_study, ".txt")
#             test_pre_dir<-paste0(project_dir, data_platform, "/", t_study, "/")
#         } else {
#             test_data_file<-paste0(project_dir, data_platform, "/", t_study, "/tg/rlog.txt")
#             test_pre_dir<-paste0(project_dir, data_platform, "/", t_study, "/tg", "/")
#         }
#         test_data<-read.table(test_data_file, sep="\t", check.names=FALSE, row.names=1)
#
#         if (study_spe == "human"){
#             ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
#             gene_name = getBM(attributes=c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"), filters="ensembl_gene_id",
#                               values=row.names(test_data), mart=ensembl_mart)
#             gene_name<-subset(gene_name, mmusculus_homolog_orthology_type=="ortholog_one2one")
#             test_data<-test_data[gene_name$ensembl_gene_id, ]
#             row.names(test_data)<-gene_name$mmusculus_homolog_ensembl_gene
#         }
#
#         # Use the markers that are both in Test and Train data
#         markers_used<-intersect(rownames(train_data_all), rownames(test_data))
#         markers_used<-intersect(rownames(marker_data), markers_used)
#
#         # Prepare the testing data
#         test_data<-test_data[markers_used, ]
#         test_data<-as.data.frame(t(test_data))
#
#         # Prepare the training data
#         train_data<-train_data_all[markers_used, ]
#         tissue_type<-str_extract(colnames(train_data), "BAT|WAT")
#         train_data<-as.data.frame(t(train_data))
#         train_data$TissueType<-as.factor(tissue_type)
#
#         # SVM
#         model_svm<-svm(TissueType ~ ., data=train_data, probability=TRUE)
#         pre_svm<-predict(model_svm, newdata=test_data, probability = TRUE)
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.svm.txt")
#         write.table(attr(pre_svm, "probabilities"), test_pre_out, sep="\t", quote=FALSE)
#
#         # naiveBayes
#         model_nB<-naiveBayes(TissueType ~ ., data=train_data, type="raw")
#         pre_nB<-predict(model_nB, newdata=test_data, type="raw")
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.nb.txt")
#         if (!identical(rownames(pre_nB), rownames(test_data))){
#             rownames(pre_nB)<-rownames(test_data)
#         }
#         write.table(pre_nB, test_pre_out, sep="\t", quote=FALSE)
#
#         # neural network
# #         train_data$TissueType<-tissue_type
# #         myform <- as.formula(paste0('TissueType ~ ', paste(names(train_data[!names(train_data) %in% 'TissueType']),
# #                                           collapse = ' + ')))
#
#         model_nnet = nnet(TissueType ~., data=train_data)
#         pre_nnet<-predict(model_nnet, newdata=test_data, type="raw")
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.nnet.txt")
#         write.table(pre_nnet, test_pre_out, sep="\t", quote=FALSE)
#
#         # random forest
#         model_rF = randomForest(TissueType ~ ., data=train_data, type="raw")
#         pre_rF<-predict(model_rF, newdata=test_data, type="prob")
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.rf.txt")
#         write.table(pre_rF, test_pre_out, sep="\t", quote=FALSE)
#
#         # recursive partioning
#         model_rpart = rpart(TissueType ~ ., data=train_data)
#         pre_rpart<-predict(model_rpart, newdata=test_data, type="prob")
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.rpart.txt")
#         write.table(pre_rF, test_pre_out, sep="\t", quote=FALSE)
#
#         # generlized linear model
#         data_x<-as.matrix(train_data[,-dim(train_data)[2]])
#         data_y<-train_data$TissueType
#
#         glm_error<-is.error(try(model_lr<-glmnet(data_x, data_y, family="multinomial", type.multinomial="grouped"), silent=TRUE))
#         #plot(lr_lamba)
#         if (glm_error){
#             pre_lr<-rep(NA, length(colnames(test_data)))
#         } else {
#             pre_res_lr<-predict(model_lr, newx=as.matrix(test_data), s=0, "response")
#         }
#         test_pre_out<-paste0(test_pre_dir, "pred.nt.SEP.glm.txt")
#         write.table(pre_res_lr[,,1], test_pre_out, sep="\t", quote=FALSE)
#     }
# }

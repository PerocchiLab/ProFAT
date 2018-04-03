# Plot the predictin accuracy.
# Last edited by Yiming on July 31 2017
# Load the data set information
library(ggplot2)
library(gplots)
library(biomaRt)
library(limma)
library(sva)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

out_dir<-paste0(project_dir, "MLResults/")

ml_method<-c("nnet")
pre_process<-NULL # NULL, scale, center. corr
if (is.null(pre_process)){
    pre_pro_txt<-"null"
} else {
    pre_pro_txt<-pre_process
}
platforms<-c("RNASeq", "MicroArray")

for (pred_method in ml_method){
    accu_file<-paste0("all.pred.wt.", pre_pro_txt, ".", pred_method, ".txt")
    for (data_platform in platforms){
        # Train and Test each study
        test_studies<-dir(paste0(project_dir, data_platform), pattern="^[mh]\\d+[MR]*", full.names=FALSE)
        for (t_study in test_studies){
            eps_file<-paste0(out_dir, t_study, ".", accu_file,".eps")
            if (file.exists(eps_file)){
               #next
            }
            if (data_platform == "RNASeq"){
                study_accu<-paste0(project_dir, "RNASeq/", t_study, "/tg/", accu_file)
                study_data<-paste0(project_dir, "RNASeq/", t_study, "/tg/", "rlog.txt")
                with_training_data<-paste0(project_dir, "RNASeq/", t_study, "/tg/", "rlog.with.training.txt")
                label_order<-paste0(project_dir, "RNASeq/", t_study, "/tg/heatmap.labels.txt")
            } else {
                study_accu<-paste0(project_dir, "MicroArray/", t_study, "/", accu_file)
                study_data<-paste0(project_dir, "MicroArray/", t_study, "/", "ND_", t_study, ".txt")
                with_training_data<-paste0(project_dir, "MicroArray/", t_study, "/", "ND_", t_study, ".with.training.txt")
                label_order<-paste0(project_dir, "MicroArray/", t_study, "/heatmap.labels.txt")
            }
            if (grepl("^m", t_study)){
                spe="mouse"
            } else if (grepl("^h", t_study)) {
                spe="human"
            }

            accu_data<-read.table(study_accu, check.names=FALSE)
            exp_data<-read.table(study_data, check.names=FALSE)
            exp_data<-get_protein_coding(exp_data, spe)
            exp_data_with_train<-read.table(with_training_data, check.names=FALSE, sep="\t")

            sample_ids<-colnames(exp_data)

            study_accu_data<-accu_data[sample_ids,]
            study_exp_data<-exp_data[, sample_ids]

            UCP1_levels<-exp_data_with_train["ENSMUSG00000031710", ]
            UCP1_levels_norm<-(UCP1_levels-min(UCP1_levels))/(max(UCP1_levels)-min(UCP1_levels))
            UCP1_levels_norm<-UCP1_levels_norm[, sample_ids]

            num_samples<-dim(study_accu_data)[1]
            lab_size<-20/num_samples
            if (lab_size>1){
                lab_size<-1
            }

            #####################################
            # Barplot figure
            #####################################
            setEPS()
            postscript(eps_file)
            par(mar=c(10,5,2,2))

            l_order<-read.table(label_order, header=TRUE, stringsAsFactors = FALSE)
            l_order<-l_order[l_order$x != "BAT" & l_order$x !="WAT",]
            study_accu_data<-study_accu_data[l_order, ]

            mp<-barplot(as.matrix(t(study_accu_data)), args.legend = list(bty="n"), beside = FALSE, ylim=c(0,1), col=c("brown", "papayawhip"),
                        las=2, cex.names=lab_size, axes = FALSE, axisnames = FALSE)
            text(mp, par("usr")[3], labels = rownames(study_accu_data), srt = 45,  adj=c(1.1, 1,1), xpd = TRUE, cex=lab_size)
            for (jj in 1:length(l_order)){
                x_axis<-1.2*jj-0.5
                y_axis<-UCP1_levels_norm[l_order[jj]]
                points(x_axis, y_axis, pch=8, cex=1.5)
            }
            dev.off()
        }
    }
}


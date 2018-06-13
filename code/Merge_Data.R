# The function is used to merge all the data from different studies.
# The intersection of all studies are merged.
library(DESeq2)
library(amap)
library(biomaRt)
library(limma)
library(sva)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

# Put all the studies data together

platforms<-c("RNASeq", "MicroArray")
#platforms<-c("MicroArray")

for (data_platform in platforms){
    ds_dir<-paste0(project_dir, data_platform)
    #studies<-dir(ds_dir, pattern="^[mh]", full.names=FALSE)
    studies<-dir(ds_dir, pattern="^[h]", full.names=FALSE)

    out_dir<-paste0(project_dir, data_platform, "/All_BeBrW")

    if (data_platform == "MicroArray"){
        study_files<-paste0(project_dir, data_platform,"/", studies, "/", "ND_", basename(studies), ".txt")
        out_f_name<-"human.All.txt"
        sum_cutoff<-0
    } else if (data_platform=="RNASeq"){
        study_files<-paste0(project_dir, data_platform,"/", studies, "/tg/rlog.txt")
        out_dir<-paste0(out_dir, "/tg")
        out_f_name<-"human.All.rlog.txt"
        sum_cutoff<-0
    }

    if (!dir.exists(out_dir)){
        dir.create(out_dir, recursive=TRUE)
    }
    out_file<-paste0(out_dir, "/", out_f_name)
    if (!file.exists(out_file)){
        gx_data<-put_together_studies(study_files, sum_cutoff)
        write.table(gx_data, out_file, sep="\t", quote=FALSE)
    } else {
        gx_data<-read.table(out_file, header=TRUE, check.names=FALSE, row.names=1)
    }
    ###################################################
    # Remove the batch effect
    ###################################################
    out_batch_file<-paste0(out_dir, "/rbatch.", out_f_name)

    #batch_VM<-gsub(".*_([MR]\\d+)\\(r\\d+\\)", "\\1", colnames(gx_data))
    batch_VM<-gsub(".*_(\\w+\\d+)$", "\\1", colnames(gx_data)) # to merge human data only. edited on 2017.09.25 by Yiming
    gx_data_batch<-ComBat(dat = gx_data, batch = batch_VM)
    write.table(gx_data_batch, out_batch_file, sep="\t", quote=FALSE)
}


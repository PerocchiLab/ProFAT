# The function is used to merge all the data from different studies.
# The intersection of all studies are merged.
library(DESeq2)
library(amap)
library(biomaRt)
library(limma)
library(sva)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

dataset_file<-paste0(project_dir, "DataInfo/RNASeq_Dataset.NoCore.txt")
data_set<-read.table(dataset_file, header=TRUE, stringsAsFactors = FALSE)
studies<-as.character(data_set$Name)

study_dir<-paste0(project_dir, "RNASeq")
studies<-dir(study_dir, pattern="^[mh]", full.names=FALSE)
studies<-c("m20R_BrW")

for (t_study in studies){
    study_file<-paste0(study_dir,"/", t_study, "/tg/seqcount.txt")
    rlog_file<-gsub("seqcount", "rlog", study_file)
    if (file.exists(rlog_file)){
        next
    }
    gx_data<-read.table(study_file, header=TRUE, check.names=FALSE, row.names=1)
    gx_data = gx_data[rowSums(gx_data)>1, ]

    sample_names<-colnames(gx_data)
    col_Data<-get_col_data(sample_names)

    dds<-DESeqDataSetFromMatrix(countData=gx_data, colData=col_Data, design=~1)
    rld <- rlog(dds, blind=TRUE)
    rlog_data<-assay(rld)

    write.table(rlog_data, gsub("seqcount", "rlog", study_file), sep="\t", quote=FALSE)

    dds <- estimateSizeFactors(dds)
    nc<-counts(dds, normalized=TRUE)
    write.table(nc, gsub("seqcount", "nc", study_file), sep="\t", quote=FALSE)
}



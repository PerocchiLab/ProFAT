# This function is used for the marker prediction.

library(caret)
library(sva)
library(stringr)
library(biomaRt)
library(limma)
library(MGFM)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

# MicroArray data
# platforms<-c("MicroArray", "RNASeq")
# # For Mouse Marker Genes.
# # Generate the significant genes for each platforms
# # Merge the significant genes across two platforms
# for (data_platform in platforms){
#     ds_dir<-paste0(project_dir, data_platform)
#
#     out_dir<-paste0(project_dir, data_platform, "/BrW_Markers")
#     if (!dir.exists(out_dir)){
#         dir.create(out_dir, recursive=TRUE)
#     }
#
#     if (data_platform == "MicroArray"){
#         # Retrieve only the datasets of BrW
#         studies<-dir(ds_dir, pattern="^[m]\\d+\\w_BrW", full.names=FALSE)
#         # Remove m09M as pvWAT_M9 is clustered together with BAT, removed for marker prediction.
#         studies<-studies[!grepl("m09M", studies)]
#
#         study_files<-paste0(project_dir, data_platform,"/", studies, "/", "ND_", basename(studies), ".txt")
#         all_data<-put_together_studies(study_files, 0)
#     } else if (data_platform == "RNASeq"){
#         studies<-dir(ds_dir, pattern="^[m]\\d+\\w_BrW", full.names=FALSE)
#         study_files<-paste0(project_dir, data_platform,"/", studies, "/tg/rlog.txt")
#         all_data<-put_together_studies(study_files, -Inf)
#     }
#
#     batch_VM<-gsub(".*_([MR]\\d+)\\(r\\d+\\)", "\\1", colnames(all_data))
#
#     tissue_type<-str_extract(colnames(all_data), "BAT|WAT")
#     tissue_type<-as.factor(tissue_type)
#
#     if (length(unique(batch_VM))>1){
#         gx_data_batch<-ComBat(dat = all_data, batch = batch_VM)
#     } else {
#         gx_data_batch<-all_data
#     }
#     write.table(all_data, paste0(out_dir, "/data_BrW.txt"), sep="\t", quote=FALSE)
#     write.table(gx_data_batch, paste0(out_dir, "/rbatch.data_BrW.txt"), sep="\t", quote=FALSE)
#
#     col_Data<-str_extract(colnames(gx_data_batch), "BAT|WAT")
#     col_Data<-as.data.frame(col_Data, row.names=colnames(gx_data_batch))
#     col_Data[,2]<-str_extract(colnames(gx_data_batch), "[MR]\\d+")
#     colnames(col_Data)<-c("Tissue", "Study")
#
#     if (length(unique(col_Data$Study))>1){
#         design<-model.matrix(~0+Tissue+Study, data=col_Data)
#     } else {
#         design<-model.matrix(~0+Tissue, data=col_Data)
#     }
#     contrast_matrix = makeContrasts(TissueBAT-TissueWAT, levels=design)
#     fit <- lmFit(all_data, design)
#     fit_cont <- contrasts.fit(fit, contrast_matrix)
#     fit_cont_eBayes <- eBayes(fit_cont)
#     num_genes = dim(fit_cont_eBayes)[1]
#     gene_list = topTable(fit_cont_eBayes, number=num_genes, sort.by="logFC", lfc=2, p.value=1e-3, confint=TRUE)
#     sig_gene_list<-get_protein_coding(gene_list, "mouse")
#
#     write.table(sig_gene_list, file=paste0(out_dir, "/sig_gene.txt"), sep="\t", quote=FALSE)
# }

# Run the marker prediction
###########################################
# Generate the combined markers for Mouse #
###########################################
# Gen_ComData(project_dir, "BrW_Markers", "mouse")
Gen_ComData(project_dir, "hBrW_Markers", "human")

# Predict the markers for each platform
# MGF_MR_Marker(project_dir)

# Predict the markers for combined data.
# MGF_MR_Marker_Com(project_dir, "BrW_Markers")
MGF_MR_Marker_Com(project_dir, "hBrW_Markers")


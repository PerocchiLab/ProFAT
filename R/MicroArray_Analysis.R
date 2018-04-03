# The function is created to retrieve the Affymetrix Microarray Data and do the Normalization
# and differential expression for Adipocyte Project.
# Authro: Yiming Cheng, Improved on 2017/01/26
# Function: This function will do the following things:
#         1. Download the data based on the GSE number and the group information
#         2. Normalize the data using quantile normalization
#         3. Use limma pakage to do the differential expression based on the group informaiton.
#         4. Save the normalized data and the differential expressed data into files.
#

#source("http://bioconductor.org/biocLite.R")
#biocLite()
#Load the necessary libraries
library(GEOquery)
#library(affy)
library(gcrma)
library(GEOmetadb)
library(limma)
library(samr)
library(LSD)
library(gplots)
library(biomaRt)
library(oligo)
library(ggplot2)

project_dir<-"/data/home/cheng/Adipocyte/"
dataset_file<-paste0(project_dir, "DataInfo/MicroArray_Dataset.NoCore.txt")

data_set<-read.table(dataset_file, header=TRUE, stringsAsFactors = FALSE)
num_dataset<-dim(data_set)[1]

for (ii in 1:num_dataset){
    dataset_info<-data_set[ii, ]

    # Load the function get_dataset_info.R
    base_dir=paste0(project_dir, "MicroArray/")  #Set working directory for download
    raw_dir = paste0(project_dir, "DataInfo/GSE/")

    ##################################################
    # Input:
    ##################################################
    GSE_Num=as.character(dataset_info$Accession)
    dataset_name = as.character(dataset_info$Name)

    samples=as.character(dataset_info$Samples)
    tissue_name=as.character(dataset_info$TissueName)
    if (grepl("^m", dataset_name)==1){
        ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl", host = "dec2015.archive.ensembl.org")
    } else if (grepl("^h", dataset_name)==1){
        ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "dec2015.archive.ensembl.org")
    }

    chip_platform = as.character(dataset_info$PlatForm)

    ##################################################
    # Retieve the data
    ##################################################
    data_dir = paste0(base_dir, dataset_name,"/")
    if (!file.exists(data_dir)){
        dir.create(data_dir, recursive=TRUE)
    }

    #Download the CEL file package for this dataset (by GSE - Geo series id)
    cel_dir = paste0(raw_dir, GSE_Num, sep = "")
    tar_file = paste0(cel_dir, "/", GSE_Num, "_RAW.tar")
    data_type = ""

    if (!file.exists(tar_file)){
        getGEOSuppFiles(GSE_Num, baseDir =raw_dir)
        # getAE("E-MTAB-758", type="raw"): this is to download the raw CEL file.
        #Unpack the CEL files
        untar(tar_file, exdir=cel_dir)
        cels = list.files(cel_dir, pattern = "cel.gz", ignore.case=TRUE)
        if (length(cels)>0){
            data_type="CEL"
            sapply(paste(cel_dir, cels, sep="/"), gunzip)
        } else { # they are not CEL file.
            try_file=paste0(cel_dir, "/", GSE_Num, "_non-normalized.txt.gz")
            if (!file.exists(try_file)){
                cat("The data are not CEL or non-normalized\n")
                return
            } else {
                data_type = "Non-Normalized"
                gunzip(try_file)
            }
        }
    }

    all_samples=strsplit(samples, ';')[[1]]
    gsm_cels = list.files(cel_dir, pattern = "CEL$", ignore.case=TRUE)

    new_ind = 1
    cels_new=c()
    for (temp_cel in all_samples){
        ind = charmatch(temp_cel, gsm_cels)
        if (is.na(ind)){
            stop("Some CEL file is missing\n")
        } else {
            cels_new[new_ind] = paste0(cel_dir, "/", gsm_cels[ind])
            new_ind = new_ind+1
        }
    }

    ##################################################
    # READ and Process the data
    ##################################################
    # Load the data and then
    # perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
    #setwd(cel_dir)

    if (grepl("_st_", chip_platform) | grepl("_hta_", chip_platform) |grepl(".db", chip_platform)){
        raw.data=read.celfiles(verbose=TRUE, filenames=cels_new) #From bioconductor
        data.rma.norm=oligo::rma(raw.data, normalize=TRUE)
    } else {
        raw.data=ReadAffy(verbose=TRUE, filenames=cels_new) #From bioconductor
        data.rma.norm=affy::rma(raw.data, normalize=TRUE)
    }

    #Get the important stuff out of the data - the expression estimates for each array
    rma_exp = exprs(data.rma.norm)

    # Name the column data
    colnames(rma_exp) = strsplit(tissue_name, ';')[[1]]

    #Format values to 5 decimal places
    rma_exp = round(rma_exp, digits=5)
    rma.avg = avereps(rma_exp) # By default, average based on row names, that is probes.
    rma.avg<-as.data.frame(rma.avg)
    rma.avg$ProbeNames<-rownames(rma.avg)

    probes=row.names(rma.avg)
    probes<-gsub("\\.\\d$", "", probes)
    row.names(rma.avg)<-probes
    if (grepl(".txt$", chip_platform)){
        probe_id = read.table(chip_platform, header=TRUE, row.names=1, sep="\t")
    } else if (grepl(".db", chip_platform)) {
        #require(chip_platform)
        library(chip_platform, character.only=TRUE)
        probe_id<-select(get(chip_platform), keys=probes, columns = c("PROBEID","ENSEMBL"))
        colnames(probe_id)<-c("ProbeNames", "ensembl_gene_id")
    } else {
        probe_id = getBM(attributes=c(chip_platform,"ensembl_gene_id"), filters=chip_platform, values=as.list(probes), mart=ensembl_mart)
        colnames(probe_id)<-c("ProbeNames", "ensembl_gene_id")
    }
    ID_data_x<-merge(rma.avg, probe_id, by.x="ProbeNames", by.y="ProbeNames")
    ID_data_x<-ID_data_x[,-1]
    ID_data_agg = aggregate(.~ensembl_gene_id, data=ID_data_x, mean)
    rownames(ID_data_agg)<-ID_data_agg$ensembl_gene_id
    ID_data_agg<-subset(ID_data_agg, select = -ensembl_gene_id)

    ID_data_x = ID_data_agg

    ##################################################
    # Write the output
    ##################################################
    #Write RMA-normalized, mapped data to file
    out_file = paste0(data_dir, "ND_", dataset_name, ".txt")
    write.table(ID_data_x, file = out_file, quote = FALSE, sep = "\t", row.names = TRUE)
}

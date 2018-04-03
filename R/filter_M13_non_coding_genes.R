# replot the PCA for microArray and RNA-Seq data
# jiang
# 2017/07/24

library(pheatmap)
library(RColorBrewer)
library(amap)
library(rafalib)
library(car)
library(ggplot2)
library(biomaRt)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
#code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
#plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))

ds_name = "mmusculus_gene_ensembl"
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=ds_name, host = "jul2015.archive.ensembl.org")


microArrayFile = list.files(path = paste0(adipocyte_dir, 'MicroArray/'), 
                            pattern = 'ND_[h|m][0-9]+M.*.txt', 
                            recursive = TRUE, full.names = TRUE)
array13 = microArrayFile[grepl('m13M', microArrayFile)] # the 13th microArray only
array13ND=  read.delim(file = array13, sep = '\t', 
                      header = TRUE, row.names = 1, check.names = FALSE)
array13ND$ensembl_gene_id = rownames(array13ND)
tmpBiotype = getBM(attributes=c("ensembl_gene_id", "gene_biotype"), filters="ensembl_gene_id",
                   values = rownames(array13ND), mart=ensembl_mart)
array13ND = merge(tmpBiotype, array13ND, by = 'ensembl_gene_id', all = FALSE)
array13ND = array13ND[array13ND$gene_biotype == 'protein_coding',]
rownames(array13ND) = array13ND$ensembl_gene_id
array13ND = array13ND[,-c(1,2)]
write.table(array13ND, file = gsub('.txt','_protein_coding.txt', array13), 
            sep = '\t', quote = FALSE)

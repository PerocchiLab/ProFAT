# Download data from GEO and EBI using esearch
# Created by Li Jiang, improved by Yiming
# 03.11.2015
library(biomaRt)
library(GEOquery)
#library(GEOmetadb)
library(stringr)
library(xlsx)
library(ArrayExpress)

source('/data/home/share/Projects/Adipocyte/rCode/Data_Download_Functions.R')
info_dir = "/data/home/share/Projects/Adipocyte/DataInfo/"

#################################################################################
# Search GEO database
#################################################################################
programmatic_data_search_dir = create_new_dir(info_dir,'programmatic_data_search/')
geo_dir = create_new_dir(programmatic_data_search_dir,'GEO/')
raw_xml_dir = create_new_dir(geo_dir,'raw_xml/')
tissue_type = c('beige','brown','white')

geo_search_keywords = c('adipocyte', paste(tissue_type, 'adipose',sep = ' AND '),
                        paste(tissue_type, 'fat',sep = ' AND '),'BAT','WAT')

## mannually download the platform information
geo_platform_dir = create_new_dir(info_dir,'geo_platform/')

## get geo platform and corresponding title
## go to the link:
#   http://www.ncbi.nlm.nih.gov/geo/browse/?view=platforms&tax=10090&zsort=date&display=20
#   http://www.ncbi.nlm.nih.gov/geo/browse/?view=platforms&tax=9606&zsort=date&display=20
#   Click "Export", then check "all search results", and "tab" delimited.
# and then click Mus musculus, then click 'all research result' choose 'tab', then save as 'mouse.tsv'
# the same method, choose 'Homo sapiens',then 'all research result', 'tab', save as 'human.tsv'

mouse_platform = read.delim(paste0(geo_platform_dir,"mouse.tsv"),header = TRUE, sep = '\t',row.names = 1)
human_platform = read.delim(paste0(geo_platform_dir,"human.tsv"),header = TRUE, sep = '\t',row.names = 1)
gse_platform = rbind(mouse_platform,human_platform)

# download the xml file by each keywords in geo_search_keywords
# download_from_geo_by_each_keyword(geo_search_keywords,raw_xml_dir)

## read all files in the raw_xml_dir and return a final data frame
summary_xml = read_all_xml_in_a_dir(raw_xml_dir)
write.table(summary_xml,file = paste0(geo_dir,'summary_xml.txt'),sep = '\t',quote = FALSE,col.names = NA)

## get biomart array name from GPL title, then can connect to ensemble id
ensembl_mart_mouse = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl",host = "jul2015.archive.ensembl.org")
array_name_mouse = listAttributes(ensembl_mart_mouse)
ensembl_mart_human = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
array_name_human = listAttributes(ensembl_mart_human)
biomart_array_name = append(array_name_human$name,array_name_mouse$name)

## get gse details and then perform further check
gse_details_dir = create_new_dir(info_dir,'geo/gse_details/')
# download all the gse details in the summary_xml
download_gse_samples_details(summary_xml$GSE,gse_details_dir)

##
# add cols to the summary_xml
for (i in 1:nrow(summary_xml)) {
    temp_gse_details = get_gse_samples_details(as.character(summary_xml$GSE[i]), gse_details_dir)
    gpl_id = as.character(summary_xml$platform_id[i])
    #summary_xml$array_name[i] =
    summary_xml$biomart[i] = from_gpl_id_to_gpl_title_then_biomart_name(gpl_id,gse_platform,biomart_array_name)
    summary_xml$contact[i] = catch_the_first_name_of_contact(temp_gse_details)
    summary_xml$total_samples[i] = nrow(temp_gse_details)
    summary_xml$total_samples_details[i] = paste_gsm_and_title_together(temp_gse_details)
    temp_gse_details_related = check_tissue_then_get_samples(temp_gse_details)
    summary_xml$related_samples[i] = nrow(temp_gse_details_related)
    summary_xml$related_samples_details[i] = paste_gsm_and_title_together(temp_gse_details_related)
    summary_xml$paired[i] = from_srp_to_paired_or_single(summary_xml$supplement[i])
    print(i)
}

final_summary_xml = summary_xml[summary_xml$related_samples > 0,]
final_summary_xml_RNAseq = final_summary_xml[grepl('SRP',final_summary_xml$supplement),]
final_summary_xml_array = final_summary_xml[!grepl('SRP',final_summary_xml$supplement),]

write.xlsx(final_summary_xml,file = paste0(geo_dir,'GEO_final_summary_xml.xlsx'),sheetName = 'all')
write.xlsx(final_summary_xml_RNAseq,file = paste0(geo_dir,'GEO_final_summary_xml.xlsx'),sheetName = 'RNAseq',append = TRUE)
write.xlsx(final_summary_xml_array,file = paste0(geo_dir,'GEO_final_summary_xml.xlsx'),sheetName = 'array',append = TRUE)

# GEO Finished
#################################################################################

#################################################################################
# extract EBI ArrayExpress and ENA
#################################################################################
ebi_dir = create_new_dir(programmatic_data_search_dir,'EBI/')
raw_txt_dir = create_new_dir(ebi_dir,'raw_txt/')

## download all the ebi datasets
download_all_ebi_from_keywords_list(geo_search_keywords,raw_txt_dir)

## read all the datasets from ebi
summary_ebi_txt = read_all_ebi_txt(raw_txt_dir)

## extract non-mouse non-human datasets and save
summary_ebi_txt_non_mouse_human = summary_ebi_txt[!(grepl('Mus|Homo',summary_ebi_txt$Organism)),]
write.table(summary_ebi_txt_non_mouse_human,file = paste0(ebi_dir,'summary_ebi_txt_non_mouse_human.txt'),
            sep = '\t',col.names = NA,quote = FALSE)

## keep only mouse and houman datasets for further analysis
summary_ebi_txt = summary_ebi_txt[grepl('Mus|Homo',summary_ebi_txt$Organism),]
write.table(summary_ebi_txt_non_mouse_human,file = paste0(ebi_dir,'summary_ebi_txt_mouse_human.txt'),
            sep = '\t',col.names = NA,quote = FALSE)

sdrf_dir = create_new_dir(info_dir,'ebi/sdrf_for_integrated_ebi_list_mRNA/')

# download all the sdrf file for the ebi
download_sdrf_file_for_all_ebi(summary_ebi_txt$Accession,sdrf_dir)

## add some cols to summary_ebi_txt
for (i in 1:nrow(summary_ebi_txt)) {
    ebi_accession = as.character(summary_ebi_txt$Accession[i])
    sdrf_file = read_single_sdrf(ebi_accession,sdrf_dir)
    sdrf_file_related = check_tissue_then_get_samples(sdrf_file)
    summary_ebi_txt$GSE[i] = match_gse(sdrf_file)
    #summary_ebi_txt$Date[i] = extract_year_from_date(summary_ebi_txt$Release.Date[i])
    summary_ebi_txt$platform[i] = check_if_ebi_array_or_seq(sdrf_file)[1]
    summary_ebi_txt$platform_name[i] = check_if_ebi_array_or_seq(sdrf_file)[2]
    summary_ebi_txt$paired[i] = paste(check_if_ebi_array_or_seq(sdrf_file)[3:length(check_if_ebi_array_or_seq(sdrf_file))],collapse = ';')
    summary_ebi_txt$total_samples[i] = nrow(sdrf_file)
    summary_ebi_txt$total_samples_details[i] = paste_df_in_one_row(get_cols_contain_factor(sdrf_file),'all')
    summary_ebi_txt$related_samples[i] = nrow(sdrf_file_related)
    summary_ebi_txt$related_samples_details[i] = paste_df_in_one_row(get_cols_contain_factor(sdrf_file_related),'related')
    print (i)
}

summary_ebi_txt = summary_ebi_txt[!grepl(';',summary_ebi_txt$GSE),]  # remove datasets with more GSE number
summary_ebi_txt = summary_ebi_txt[!grepl(';',summary_ebi_txt$platform_name),] # remove datasets with more platform
summary_ebi_txt = summary_ebi_txt[!grepl(';',summary_ebi_txt$paired),]
dim(summary_ebi_txt) # 420*14
### check if GEOD in GEO exported list
summary_ebi_txt_uniq_in_ebi = remove_ebi_dataset_which_can_be_found_in_geo_search(summary_ebi_txt,summary_xml)

dim(summary_ebi_txt_uniq_in_ebi) # 216*14

adf_dir = create_new_dir(info_dir,'ebi/adf_for_selected_ebi_datasets/')
for (i in 1:nrow(summary_ebi_txt_uniq_in_ebi)) {
    platform_name = as.character(summary_ebi_txt_uniq_in_ebi$platform_name[i])
    accession_id = as.character(summary_ebi_txt_uniq_in_ebi$Accession[i])
    summary_ebi_txt_uniq_in_ebi$biomart[i] = from_ebi_platform_name_to_biomart(platform_name,
                                                                               accession_id,gse_platform,
                                                                               biomart_array_name,adf_dir)
    summary_ebi_txt_uniq_in_ebi$Pubmed[i] = queryAE_info_pubmed(accession_id)
    print (i)
}

# remove datasets that from non Affymetrix or HiSeq platform
summary_ebi_txt_uniq_in_ebi_affy_hiseq = summary_ebi_txt_uniq_in_ebi[(!grepl('NoAffy',summary_ebi_txt_uniq_in_ebi$biomart)),]
dim(summary_ebi_txt_uniq_in_ebi_affy_hiseq) # 30*16
summary_ebi_txt_uniq_in_ebi_affy_hiseq = summary_ebi_txt_uniq_in_ebi_affy_hiseq[(
    !grepl('Not|need_test|NoArrayName|Name',summary_ebi_txt_uniq_in_ebi_affy_hiseq$biomart)),]

## I have checked the remove dataset with biomart information 'Not found','need_test','NoArrayName'
# E-GEOD-50574 is a miRNA array

summary_ebi_txt_uniq_in_ebi_affy_hiseq =
    summary_ebi_txt_uniq_in_ebi_affy_hiseq[summary_ebi_txt_uniq_in_ebi_affy_hiseq$related_samples >0,]

# remove dataset without related samples
summary_ebi_txt_uniq_in_ebi_affy_hiseq =
    summary_ebi_txt_uniq_in_ebi_affy_hiseq[(!grepl('repeated',summary_ebi_txt_uniq_in_ebi_affy_hiseq$total_samples_details)),]

#### combine two datasets:
summary_xml_selected_cols = final_summary_xml[,c(
    "GSE","supplement",
    "title","species","pubmed",
    "date","biomart","platform_id",
    "paired",
    "total_samples","total_samples_details",
    "related_samples","related_samples_details"
)]

# change the second cols the col name is 'data_accession'
colnames(summary_xml_selected_cols)[2] = 'data_accession'
for (i in 1:nrow(summary_xml_selected_cols)) {
    gse = as.character(summary_xml_selected_cols$GSE[i])
    temp_data_accession = as.character(summary_xml_selected_cols$data_accession[i])
    if (grepl('SRP',temp_data_accession)) {
    }else{
        summary_xml_selected_cols$data_accession[i] = gse
    }
    print (i)
}

## ebi output   # summary_ebi_txt_uniq_in_ebi_affy_hiseq
# add another col which is same to the Accession col
summary_ebi_txt_uniq_in_ebi_affy_hiseq$data_accession = summary_ebi_txt_uniq_in_ebi_affy_hiseq$Accession
ebi_summary_selected_cols = summary_ebi_txt_uniq_in_ebi_affy_hiseq[,c(
    "Accession","data_accession",
    "Title","Organism","Pubmed",
    "Release.Date","biomart","platform_name",
    "paired",
    "total_samples","total_samples_details",
    "related_samples","related_samples_details"
)]
# make the two datasets with the same col names
colnames(summary_xml_selected_cols) = colnames(ebi_summary_selected_cols)
all_geo_and_ebi = rbind(summary_xml_selected_cols,ebi_summary_selected_cols)
all_geo_and_ebi = all_geo_and_ebi[!grepl(';',all_geo_and_ebi$paired),]

for (i in 1:nrow(all_geo_and_ebi)) {
    paired = all_geo_and_ebi$paired[i]
    if (grepl('[aA]rray',paired)) {
        all_geo_and_ebi$paired[i] = -1
        #seq_type = append(seq_type,'M')
    }else{
        #seq_type = append(seq_type,'R')
    }
}

colnames(all_geo_and_ebi)[grep('pair',colnames(all_geo_and_ebi))] = "sequence_type"

all_geo_and_ebi_mouse = all_geo_and_ebi[grepl('Mus',all_geo_and_ebi$Organism),]
all_geo_and_ebi_human = all_geo_and_ebi[grepl('Homo',all_geo_and_ebi$Organism),]
rownames(all_geo_and_ebi_mouse) = paste0('m',formatC(1:nrow(all_geo_and_ebi_mouse), width = 3, flag = "0"),
                                         get_seq_type_list(all_geo_and_ebi_mouse))
rownames(all_geo_and_ebi_human) = paste0('h',formatC(1:nrow(all_geo_and_ebi_human), width = 3, flag = "0"),
                                         get_seq_type_list(all_geo_and_ebi_human))

date_str<-format(Sys.Date(), format="%Y.%m.%d")
out_file_name<-paste0(programmatic_data_search_dir,'all_geo_ebi_', date_str, '.xlsx')

write.xlsx(all_geo_and_ebi_mouse,file = out_file_name,sheetName = 'Mouse')
write.xlsx(all_geo_and_ebi_human,file = out_file_name,sheetName = 'Human',append = TRUE)

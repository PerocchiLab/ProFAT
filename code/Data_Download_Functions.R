from_string_to_valid_percentage = function(a_line_as_string){
    split_list = unlist(strsplit(a_line_as_string,'\t|;'))
    cluster_number_list = as.integer(split_list[2:4])
    valid_percentage = cluster_number_list[2]/cluster_number_list[1]
    #correct_percentage = cluster_number_list[3]/cluster_number_list[2]
    return (valid_percentage)
}

from_string_to_correct_percentage = function(a_line_as_string){
    split_list = unlist(strsplit(a_line_as_string,'\t|;'))
    cluster_number_list = as.integer(split_list[2:4])
    #valid_percentage = cluster_number_list[2]/cluster_number_list[1]
    correct_percentage = cluster_number_list[3]/cluster_number_list[2]
    #return (valid_percentage)
    return(correct_percentage)
}

from_string_to_row_name = function(a_line_as_string){
    split_list = unlist(strsplit(a_line_as_string,'\t|;'))
    row_name = split_list[1]
    return (row_name)
}


catch_percentage_row = function(keyword_for_search,cluster_line){
    line_number = grep(keyword_for_search,cluster_line)
    if (nchar(cluster_line[line_number + 4])>0){
        temp_line = cluster_line[(line_number + 1):(line_number + 4)]
    }else
    { temp_line = cluster_line[(line_number + 1):(line_number + 3)]}
    return(temp_line)
}


catch_percentage_row_number = function(keyword_for_search,cluster_line){
    line_number = grep(keyword_for_search,cluster_line)
    if (nchar(cluster_line[line_number + 4])>0){
        #temp_line = cluster_line[(line_number + 1):(line_number + 4)]
        row_number = 4
    }else
    {
        #temp_line = cluster_line[(line_number + 1):(line_number + 3)]
        row_number =3 }
    return(row_number)
}

bar_plot_valid = function(percentage_df){
    RR = 'lightblue'
    MM = 'green'
    All = "blue"
    MR = 'red'
    mybarcol = 'gray20'
    plot_color = character()
    for (i in 1:nrow(percentage_df)){
        temp_color = get(rownames(percentage_df)[i])
        plot_color = append(plot_color,temp_color)
    }
    barplot(as.matrix(percentage_df), beside = TRUE,col = plot_color,
            legend = rownames(percentage_df), ylim = c(0,1.0), xlim = c(0,5.8*ncol(percentage_df)),
            main = "Valid Percentage", font.main = 4, col.sub = mybarcol,
            cex.names = 1.5,args.legend = list(bty="n"))
}

bar_plot_correct = function(percentage_df){
    RR = 'lightblue'
    MM = 'green'
    All = "blue"
    MR = 'red'
    mybarcol = 'gray20'
    plot_color = character()
    for (i in 1:nrow(percentage_df)){
        temp_color = get(rownames(percentage_df)[i])
        plot_color = append(plot_color,temp_color)
    }
    barplot(as.matrix(percentage_df), beside = TRUE,col = plot_color,
            legend = rownames(percentage_df), ylim = c(0,1.0), xlim = c(0,5.8*ncol(percentage_df)),
            main = "Correct Percentage", font.main = 4, col.sub = mybarcol,
            cex.names = 1.5,args.legend = list(bty="n"))
}

bar_plot = function(percentage_df){
    RR = 'lightblue'
    MM = 'green'
    All = "blue"
    MR = 'red'
    mybarcol = 'gray20'
    plot_color = character()
    for (i in 1:nrow(percentage_df)){
        temp_color = get(rownames(percentage_df)[i])
        plot_color = append(plot_color,temp_color)
    }
    barplot(as.matrix(percentage_df), beside = TRUE,col = plot_color,
            legend = rownames(percentage_df), ylim = c(0,1.0), xlim = c(0,5.8*ncol(percentage_df)),
            main = "All Parameters", font.main = 4, col.sub = mybarcol,
            cex.names = 1,args.legend = list(bty="n"))
}

create_new_dir = function(dir_path,dir_name) {
    # give the dir path to a variable, if the path does not exist, create it
    new_dir = paste0(dir_path,dir_name)
    if (!file.exists(new_dir)) {
        dir.create(new_dir, recursive = TRUE)
    }
    return (new_dir)
}

turn_keywords_to_file_name = function(the_keyword) {
    # give a string contain space or not,then connect the string with _
    new_list = unlist(strsplit(the_keyword,'\\s'))
    file_name = paste(new_list,collapse = '_')
    return(file_name)
}

download_from_geo_by_each_keyword = function(keyword_list,file_to_save) {
    # download data by esearch and esummary
    direct_dir = '/data/home/share/Tools/edirect/'
    for (i in 1:length(keyword_list)) {
        keyword = keyword_list[i]
        # new_command = paste0(direct_dir,"esearch -db gds -query '",keyword,
        #                     "' | ",direct_dir,"esummary -format -mode xml -db gds >",file_to_save,
        #                    turn_keywords_to_file_name(keyword),
        #                   '.txt')
        # with '-format' parameters in esummary, the result will be text format rather than xml format
        # we use xml format because xml format contain more information than the txt format
        new_command = paste0(
            direct_dir,"esearch -db gds -query '",keyword,
            "' | ",direct_dir,"esummary -mode xml -db gds > ",file_to_save,
            turn_keywords_to_file_name(keyword),
            '.xml'
        )
        system(new_command)
    }
}


extract_gds_from_lines = function(lines_file) {
    if (sum(grepl('GDS\\d+',lines_file))) {
        temp_lines = lines_file[grepl('GDS\\d+',lines_file)]
        gds = unlist(str_split(temp_lines[1],'>|<'))[3]
    }else{
        gds = 'GDS not found'
    }
    return (gds)
}

extract_gse_from_lines = function(lines_file) {
    if (sum(grepl('<GSE>',lines_file))) {
        temp_lines = lines_file[grepl('<GSE>',lines_file)]
        gse = paste(unlist(str_split(temp_lines[1],'>|<'))[2:3],collapse = "")
    }else{
        gse = 'GSE not found'
    }
    return (gse)
}

extract_gpl_from_lines = function(lines_file) {
    if (sum(grepl('<GPL>',lines_file))) {
        temp_lines = lines_file[grepl('<GPL>',lines_file)]
        gpl = paste(unlist(str_split(temp_lines[1],'>|<'))[2:3],collapse = "")
    }else{
        gpl = 'GPL not found'
    }
    return (gpl)
}

extract_gse_title_from_lines = function(lines_file) {
    if (sum(grepl('<title>',lines_file))) {
        temp_lines = lines_file[grepl('<title>',lines_file)]
        gse_title = unlist(str_split(temp_lines[1],'>|<'))[3]
    }else{
        gse_title = 'title not found'
    }
    return (gse_title)
}


extract_species_from_lines = function(lines_file) {
    if (sum(grepl('<taxon>',lines_file))) {
        temp_lines = lines_file[grepl('<taxon>',lines_file)]
        species = unlist(str_split(temp_lines[1],'>|<'))[3]
    }else{
        species = 'taxon not found'
    }
    return (species)
}



extract_pubmedid_from_lines = function(lines_file) {
    if (sum(grepl('PubMedIds',lines_file))) {
        index_pubmedids = grep('PubMedIds',lines_file)
        temp_lines = lines_file[index_pubmedids[1] + 1]
        pubmedid = as.integer(unlist(str_split(temp_lines[1],'>|<'))[3])
    }else{
        pubmedid = 'PumMedIds not found'
    }
    return (pubmedid)
}



extract_date_from_lines = function(lines_file) {
    if (sum(grepl('PDAT',lines_file))) {
        temp_lines = lines_file[grepl('PDAT',lines_file)]
        date = unlist(str_split(temp_lines[1],'>|<'))[3]
        date = unlist(strsplit(date,'/'))[1]
    }else{
        date = 'PDAT not found'
    }
    return (date)
}



extract_platform_title_from_lines = function(lines_file) {
    if (sum(grepl('PlatformTitle',lines_file))) {
        temp_lines = lines_file[grepl('PlatformTitle',lines_file)]
        platform_title = unlist(str_split(temp_lines[1],'\\[|\\]'))[2]
        platform_title = gsub('-','_',platform_title)
        platform_title = tolower(platform_title)
    }else{
        platform_title = 'PlatformTitle not found'
    }
    return (platform_title)
}

check_if_data_is_gsm_gpl = function(lines_file) {
    if (sum(grepl('Accession',lines_file))) {
        temp_lines = lines_file[grepl('Accession',lines_file)]
        if (sum(grepl('GSM|GPL',temp_lines[1]))) {
            gsm = "YES"
        }else{
            gsm = "NO"
        }
    }else{
        gsm = 'Not found'
    }
    return (gsm)
}


extract_supplement_type_from_lines = function(lines_file) {
    if (sum(grepl('SRP',lines_file))) {
        temp_lines = lines_file[grepl('SRP',lines_file)]
        supplement = unlist(str_split(temp_lines[1],'>|<'))[3]

    }else{
        if (sum(grepl('suppFile',lines_file))) {
            temp_lines = lines_file[grepl('suppFile',lines_file)]
            supplement = unlist(str_split(temp_lines[1],'>|<'))[3]
        }else{
            supplement = 'suppFile not found'
        }
    }

    return (supplement)
}




check_if_gds_type_is_expression = function(lines_file) {
    if (sum(grepl('gdsType',lines_file))) {
        temp_lines = lines_file[grepl('gdsType',lines_file)]
        if (sum(grepl('Expression',temp_lines[1]))) {
            expression = "YES"
        }else{
            expression = "NO"
        }
    }else{
        expression = 'Not found'
    }
    return (expression)
}




## read xml file
read_xml_file = function(xml_file_dir,xml_file_name) {
    temp_xml = paste0(xml_file_dir,xml_file_name)
    conn = file(temp_xml,open = "r")
    linn = readLines(conn) # read xml file line by line
    #index_DocumentSummary = grep('DocumentSummary>',linn)
    index_DocumentSummary_start = grep('<DocumentSummary>',linn)
    index_DocumentSummary_end = grep('</DocumentSummary>',linn)
    #index_DocumentSummary = index_DocumentSummary[2:length(index_DocumentSummary)-1]
    # remove the 1st and last match, because they are not for single datasets,the are title for all the file,
    #the mathch DocumentSummarySet, but using 'DocumentSummary>' we can only match the wanted part
    df_for_all_gds = as.data.frame(matrix(NA,0,9))
    colnames(df_for_all_gds) = c(
        'GDS','GSE','title','species',
        'pubmed','date','platform_id','platform','supplement'
    )
    for (i in 1:length(index_DocumentSummary_start)) {
        temp_linn = linn[index_DocumentSummary_start[i]:index_DocumentSummary_end[i]]
        if (check_if_data_is_gsm_gpl(temp_linn) != "YES") {
            if (sum(grepl(
                'CEL|SRP|BED',extract_supplement_type_from_lines(temp_linn)
            ))) {
                if (check_if_gds_type_is_expression(temp_linn) != 'NO') {
                    if (grepl('Mus|Homo',extract_species_from_lines(temp_linn))) {
                        if (sum(!grepl(
                            ';',extract_species_from_lines(temp_linn)
                        ))) {
                            if (!(grepl(
                                ';',extract_gpl_from_lines(temp_linn)
                            ))) {
                                df_for_all_gds[i,1] = extract_gds_from_lines(temp_linn)
                                df_for_all_gds[i,2] = extract_gse_from_lines(temp_linn)
                                df_for_all_gds[i,3] = extract_gse_title_from_lines(temp_linn)
                                df_for_all_gds[i,4] = extract_species_from_lines(temp_linn)
                                df_for_all_gds[i,5] = extract_pubmedid_from_lines(temp_linn)
                                df_for_all_gds[i,6] = extract_date_from_lines(temp_linn)
                                df_for_all_gds[i,7] = extract_gpl_from_lines(temp_linn)
                                df_for_all_gds[i,8] = extract_platform_title_from_lines(temp_linn)
                                df_for_all_gds[i,9] = extract_supplement_type_from_lines(temp_linn)
                            }
                        }
                    }
                }
            }
        }
        #print (i)
    }
    df_for_all_gds = df_for_all_gds[!duplicated(df_for_all_gds),]
    df_for_all_gds = df_for_all_gds[as.logical(rowSums(!is.na(df_for_all_gds))),]
    close(conn)
    return (df_for_all_gds)
}



read_all_xml_in_a_dir = function(dir_name) {
    all_files = list.files(dir_name,pattern = '.xml')
    final_df = as.data.frame(matrix(NA,0,9))
    colnames(final_df) = c(
        'GDS','GSE','title','species',
        'pubmed','date','platform_id','platform','supplement'
    )
    for (i in 1:length(all_files)) {
        temp_file = read_xml_file(dir_name,all_files[i])
        final_df = rbind(final_df,temp_file)
        print (i)
    }
    final_df = final_df[!duplicated(final_df),]
    final_df = final_df[!duplicated(final_df$GSE),]
    return (final_df)
}




# FUNCTION : catch the array name for a given GPL title
# title example:"[HG-U219_Hs_ENTREZG_18] Affymetrix Human Genome U219 Array [CDF: Brainarray HGU219_Hs_ENTREZG_v18]"
catch_array_name = function(the_string) {
    the_string = as.character(the_string)
    if (grepl('\\[',the_string)) {
        temp_name = unlist(str_split(the_string, '\\]|\\[', n = Inf))[2]
        name = gsub('-','_',temp_name)
        name = tolower(name)
        return(name)
    }
    else{
        if (grepl('HiSeq',the_string)) {
            name = the_string
        }else{
            name = 'No gpl Name'
        }
        return(name)
    }
}


# from platform title name to biomart array name
from_platform_name_to_biomart_name = function(the_string,biomart_array_character) {
    gse_array_name = catch_array_name(the_string)
    if (gse_array_name != 'No gpl Name') {
        if (sum(grepl('HiSeq',gse_array_name))) {
            return (gse_array_name)
        } else{
            if (sum(grepl(gse_array_name,biomart_array_character))) {
                biomart_name = paste(biomart_array_character[grep(gse_array_name,biomart_array_character)],collapse = ';')
                return(biomart_name)
            }else {
                return("No biomart Name")
            }
        }
    }else {
        return('No gpl Name')
    }
}

## function

from_gpl_id_to_gpl_title_then_biomart_name = function(single_gpl_id,gpl_df,biomart_array_name_list) {
    single_gpl_id = as.character(single_gpl_id)
    gpl_title = as.character(gpl_df[single_gpl_id,1])
    if (sum(grepl('[aA]ffy|HiSeq|Seq',gpl_title)) > 0) {
        biomart_name = from_platform_name_to_biomart_name(gpl_title,biomart_array_name_list)
    }else{
        biomart_name = 'NoAffyOrSeq'
    }
    return(biomart_name)
}


## download all the gse detail for a give list of gse number,then save then in a  given dir
download_gse_samples_details = function(gse_number_list, gse_details_save_dir) {
    gse_number_list = as.character(gse_number_list)
    for (i in 1:length(gse_number_list)) {
        gse_file = paste0(gse_details_save_dir, gse_number_list[i], ".txt")
        if (file.exists(gse_file)) {
            #tmp_info = read.delim(gse_file,header = TRUE,sep = '\t')
        } else {
            tmp_info = pData(getGEO(gse_number_list[i], GSEMatrix=FALSE, getGPL=FALSE)[[1]])
            tmp_info$GSE_accession = gse_number_list[i]
            write.table(
                tmp_info,file = paste0(gse_details_save_dir, gse_number_list[i], ".txt"),
                sep = '\t',quote = FALSE,col.names = NA
            )
        }
    }
}

# read gse samples detail,if file does not exist, download and save
get_gse_samples_details = function(gse_number, gse_details_save_dir) {
    gse_file = paste0(gse_details_save_dir, gse_number, ".txt")
    if (file.exists(gse_file)) {
        try_tmp_info = try(read.delim(gse_file,header = TRUE,sep = '\t'))
        if (sum(grepl('Error',try_tmp_info))) {
            tmp_info = pData(getGEO(gse_number)[[1]])
            tmp_info$GSE_accession = gse_number
        }else{
            tmp_info = read.delim(gse_file,header = TRUE,sep = '\t')
        }

    } else {
        tmp_info = pData(getGEO(gse_number)[[1]])
        tmp_info$GSE_accession = gse_number
        write.table(
            tmp_info,file = paste0(gse_details_save_dir, gse_number, ".txt"),
            sep = '\t',quote = FALSE,col.names = NA
        )
    }
    return(tmp_info)
}


### catch first name of the contact, from the temp_gse_details file
##### FUNCTION

catch_the_first_name_of_contact = function(gse_details_df) {
    if (sum(grepl('contact_name',colnames(gse_details_df)))) {
        the_string = as.character(levels(gse_details_df$contact_name))
        first_name = unlist(str_split(the_string, '\\s|,', n = Inf))[1]
    }else{
        first_name = 'No contact_name'
    }
    return(first_name)
}


paste_gsm_and_title_together = function(gse_details_df) {
    if (nrow(gse_details_df) > 50) {
        gsm_title = 'samples > 50'
    }else{
        if (nrow(gse_details_df) > 0) {
            gsm_title = character()
            for (i in 1:nrow(gse_details_df)) {
                temp_gsm_title = paste0(gse_details_df$geo_accession[i],'\t',gse_details_df$title[i],'\n')
                gsm_title = paste0(gsm_title,temp_gsm_title)
            }
        }else{
            gsm_title = '0 sample'
        }
    }
    return(gsm_title)
}


#########
## Function: check tissues then return related samples in a dataframe, rather than the number of samples
#####################################
check_BrBeW_samples = function(df) {
    samples_having_BrBeW = as.logical(rowSums(apply(df, 2, function(x)
        grepl('[wW]hite|[Bb]rown|[Bb]eige|BA|WA', x))))
    new_df_contain_tissue = df[samples_having_BrBeW,]
    return(new_df_contain_tissue)
    # this will return rows that contain theses key words
}

## Function : give an data frame in human
check_adipose_samples = function(df) {
    samples_having_adipose = as.logical(rowSums(apply(df, 2, function(x)
        grepl('[Aa]dipose|[Aa]dipocyte|[Ff]at|BA|WA', x))))
    new_df_contain_tissue = df[samples_having_adipose,]
    return(new_df_contain_tissue)
}

####
check_tissue_then_get_samples = function(df) {
    if (nrow(df) > 1) {
        mouse_factor = as.logical(rowSums(apply(df,2,function(x) {
            grepl('Mus\\smusculus', x)
        }))) # which samples are mouses
        if (sum(mouse_factor)) {
            new_df_mouse = check_BrBeW_samples(df[mouse_factor,])
        }
        else {
            new_df_mouse = as.data.frame(matrix(NA,0,ncol(df)))
        }

        if (sum(!(mouse_factor))) {
            new_df_non_mouse = check_adipose_samples(df[!(mouse_factor),])
        }
        else{
            new_df_non_mouse = as.data.frame(matrix(NA,0,ncol(df)))
        }
        related_samples = rbind(new_df_mouse,new_df_non_mouse)
    }else{
        related_samples = as.data.frame(matrix(NA,0,ncol(df)))
    }

    return(related_samples)
}

from_srp_to_paired_or_single = function(SRP_id) {
    if (sum(grepl('SRP',SRP_id))) {
        tool_dir = '/data/home/share/Tools/edirect/'
        new_command = paste0(
            tool_dir,"esearch -db sra -query '", SRP_id,"' | ",tool_dir,"efetch -format docsum"
        )
        screen_output = system(new_command,intern = TRUE)
        layout_info = unique(str_extract(screen_output,'<LIBRARY_LAYOUT>\\s+<(\\w+)/'))
        if (sum(grepl('PAIRED',layout_info))) {
            if (sum(grepl('SINGLE',layout_info))) {
                paired = c(0,1)
            }else{
                paired = 1
            }
        }else{
            if (sum(grepl('SINGLE',layout_info))) {
                paired = 0
            }else{
                paired = -1
            }
        }
        return (paste0(paired,collapse = ';'))
    }else{
        return ('Array')
    }
}

####################################################################3
# EBI
######################################################################
turn_keywords_to_ebi_search_format = function(the_keyword) {
    # give a string contain space or not,then connect the string with _
    new_list = unlist(strsplit(the_keyword,'\\s'))
    file_name = paste(new_list,collapse = '+')
    return(file_name)
}


download_single_ebi = function(the_keyword, dir_to_save) {
    link = paste0(
        'https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?query=',
        turn_keywords_to_ebi_search_format(the_keyword)
    )
    wget_command = paste0('wget -O ',dir_to_save,turn_keywords_to_file_name(the_keyword),'.txt ',link)
    system(wget_command)
}

## download dataset with all the keywords
download_all_ebi_from_keywords_list = function(keywords_list,dir_to_save) {
    for (i in 1:length(keywords_list)) {
        download_single_ebi(keywords_list[i],dir_to_save)
    }
}

## chcek if datasets are for message RNA
check_if_mRNA_data = function(df) {
    platform_factor = grepl('RNA-seq\\sof\\scoding\\sRNA|transcription', df$Type)
    temp_df = df[platform_factor,]
    return(temp_df)
}


## read one txt file
read_single_ebi_txt = function(file_dir,file_name) {
    single_txt = read.delim(paste0(file_dir,file_name), header = TRUE, sep =
                                "\t")
    single_txt = single_txt[!grepl(',',single_txt$Organism),] # keep datasets than only have one organism,
    # different organisms are seperated by ','
    single_txt = single_txt[single_txt$Assays > 1,]
    single_txt = check_if_mRNA_data(single_txt)
    single_txt = single_txt[,1:6]
    return(single_txt)
}

## read all .txt files in a dir

read_all_ebi_txt = function(file_dir) {
    files_list = list.files(file_dir,pattern = '.txt')
    for (i in 1:length(files_list)) {
        if (!exists('all_ebi_txt')) {
            all_ebi_txt = read_single_ebi_txt(file_dir,files_list[i])
        }else{
            temp_single_txt = read_single_ebi_txt(file_dir,files_list[i])
            all_ebi_txt = rbind(all_ebi_txt,temp_single_txt)
        }
    }
    all_ebi_txt = all_ebi_txt[!duplicated(all_ebi_txt),]
    all_ebi_txt = all_ebi_txt[!(duplicated(all_ebi_txt$Accession)),]
    return(all_ebi_txt)
}



match_gse = function(sdrf_df) {
    temp_gse_list = apply(sdrf_df, 2, function(x) {
        if (sum(!is.na(x))) {
            ifelse(gregexpr("GSE\\d+", as.character(x)),return(unlist(unique(
                regmatches(x, gregexpr("GSE\\d+", as.character(x)))
            ))),)
        }
    })
    temp_gse = unique(unlist(temp_gse_list))
    if (!(length(temp_gse) == 0))  {
        temp_gse_list = paste(temp_gse,collapse = ';')
        return(temp_gse_list)
    }
    else {
        return ('NO')
    }
}

download_sdrf_file_for_all_ebi = function(ebi_accession_list,sdrf_save_dir) {
    ebi_accession_list = as.character(ebi_accession_list)
    for (i in 1:length(ebi_accession_list)) {
        sdrf_file = paste0(sdrf_save_dir,ebi_accession_list[i],'.sdrf.txt')
        if (!file.exists(sdrf_file)) {
            system(
                paste0(
                    'wget https://www.ebi.ac.uk/arrayexpress/files/',ebi_accession_list[i],'/',
                    ebi_accession_list[i],'.sdrf.txt ',
                    '-O ',sdrf_save_dir,ebi_accession_list[i],'.sdrf.txt'
                )
            )
        }
        print (i)
    }
}

# read_single sdrf file for each ebi accession

read_single_sdrf = function(single_ebi_accession,sdrf_save_dir) {
    single_ebi_accession = as.character(single_ebi_accession)
    temp = try(read.delim(
        paste0(sdrf_save_dir,single_ebi_accession,'.sdrf.txt'),sep = '\t',header = TRUE
    ))
    if (inherits(temp, 'try-error')) {
        system(
            paste0(
                'wget https://www.ebi.ac.uk/arrayexpress/files/',single_ebi_accession,'/',
                single_ebi_accession,'.hyb.sdrf.txt ',
                '-O ',sdrf_save_dir,single_ebi_accession,'.sdrf.txt'
            )
        )
    }
    temp_sdrf = read.delim(
        paste0(sdrf_save_dir,single_ebi_accession,'.sdrf.txt'),sep = '\t',header = TRUE
    )
    return(temp_sdrf)
}

extract_year_from_date = function(the_date) {
    the_date = as.character(the_date)
    date_list = unlist(strsplit(the_date,'-|/'))
    year = date_list[1]
    return(year)
}

check_if_ebi_array_or_seq = function(sdrf_temp_file) {
    if (sum(grepl('Array.Design.REF',colnames(sdrf_temp_file)))) {
        data_type = 'microArray'
        type_name = paste(as.character(unique(
            sdrf_temp_file$Array.Design.REF
        )),collapse = ';')
        paired = 'array'
        data_type_and_ebi_platform_name = c(data_type,type_name,paired)
        return(data_type_and_ebi_platform_name)
    }
    else{
        if (sum(grepl('INSTRUMENT_MODEL',colnames(sdrf_temp_file)))) {
            data_type = 'RNA-Seq'
            type_name = paste(as.character(unique(sdrf_temp_file[,colnames(sdrf_temp_file)[grepl('INSTRUMENT_MODEL',colnames(sdrf_temp_file))]])),collapse = ';')
            if (sum(grepl('LIBRARY_LAYOUT',colnames(sdrf_temp_file)))) {
                paired_name = paste(as.character(unique(sdrf_temp_file[,colnames(sdrf_temp_file)[grepl('LIBRARY_LAYOUT',colnames(sdrf_temp_file))]])),collapse = ';')
                if (sum(grepl('pair',paired_name,ignore.case = TRUE))) {
                    if (sum(grepl('single',paired_name,ignore.case = TRUE))) {
                        paired = c(0,1)
                    }else{
                        paired = 1
                    }
                }else{
                    if (sum(grepl('single',paired_name,ignore.case = TRUE))) {
                        paired = 0
                    }else{
                        paired = -1
                    }
                }
            }else{
                paired = -1
            }
            data_type_and_ebi_platform_name = c(data_type,type_name,paired)
            return(data_type_and_ebi_platform_name)
        }else{
            if (sum(grepl('LIBRARY_STRATEGY',colnames(sdrf_temp_file)))) {
                data_type = 'RNA-Seq'
                type_name = paste(as.character(unique(sdrf_temp_file[,colnames(sdrf_temp_file)[grepl('LIBRARY_STRATEGY',colnames(sdrf_temp_file))]])),
                                  collapse = ';')
                if (sum(grepl('LIBRARY_LAYOUT',colnames(sdrf_temp_file)))) {
                    paired_name = paste(as.character(unique(sdrf_temp_file[,colnames(sdrf_temp_file)[grepl('LIBRARY_LAYOUT',colnames(sdrf_temp_file))]])),collapse = ';')
                    if (sum(grepl('pair',paired_name,ignore.case = TRUE))) {
                        if (sum(grepl('single',paired_name,ignore.case = TRUE))) {
                            paired = c(0,1)
                        }else{
                            paired = 1
                        }
                    }else{
                        if (sum(grepl('single',paired_name,ignore.case = TRUE))) {
                            paired = 0
                        }else{
                            paired = -1
                        }
                    }
                }
                if (type_name == 'RNA-Seq') {
                    if (sum(grepl(
                        'Platform_title',colnames(sdrf_temp_file)
                    ))) {
                        type_name = paste(as.character(unique(
                            sdrf_temp_file[,colnames(sdrf_temp_file)[grep('Platform_title',colnames(sdrf_temp_file))]]
                        )),
                        collapse = ';')
                    }
                }
                data_type_and_ebi_platform_name = c(data_type,type_name,paired)
                return(data_type_and_ebi_platform_name)

            }else {
                data_type = 'need_test'
                type_name = 'need_test'
                paired = 'need_test'
                data_type_and_ebi_platform_name = c(data_type,type_name,paired)
                return(data_type_and_ebi_platform_name)
            }

        }
    }
}


catch_the_first_line_adf_file = function(file_name) {
    temp_adf = file(file_name,open = 'r')
    temp_adf_all_line = readLines(temp_adf)
    the_1st_string = temp_adf_all_line[1]
    return(the_1st_string)
}

## given an line or string, catch the conten into [] and then convert into lower
catch_array_name = function(the_string) {
    the_string = as.character(the_string)
    if (grepl('\\[',the_string)) {
        temp_name = unlist(str_split(the_string, '\\]|\\[', n = Inf))[2]
        name = gsub('-','_',temp_name)
        name = tolower(name)
        return(name)
    }
    else{
        return('RNA-Seq')
    }
}


from_array_name_to_biomart_name = function(adf_file,biomart_array_name) {
    gse_array_name = catch_array_name(catch_the_first_line_adf_file(adf_file))
    if (!gse_array_name == 'RNA-Seq')
    {
        if (sum(grepl(gse_array_name,biomart_array_name))) {
            biomart_name = paste(biomart_array_name[grep(gse_array_name,biomart_array_name)],collapse = ';')
            return(biomart_name)
        }
        else {
            return("Not Found")
        }
    }
    else {
        return('RNA-Seq')
    }
}

from_ebi_GEOD_to_gpl = function(ebi_GEOD) {
    ebi_GEOD = as.character(ebi_GEOD)
    geod_list = unlist(strsplit(ebi_GEOD,'-'))
    gpl_name = paste0('GPL',geod_list[3])
    return(gpl_name)
}






catch_the_first_line_adf_file = function(file_name) {
    temp_adf = file(file_name,open = 'r')
    temp_adf_all_line = readLines(temp_adf)
    the_1st_string = temp_adf_all_line[1]
    close(temp_adf)
    return(the_1st_string)
}

catch_array_name = function(the_string) {
    the_string = as.character(the_string)
    if (sum(grepl('Affymetrix',the_string,ignore.case = TRUE))) {
        if (grepl('\\[',the_string)) {
            temp_name = unlist(str_split(the_string, '\\]|\\[', n = Inf))[2]
            name = gsub('-','_',temp_name)
            name = tolower(name)
            return(name)
        }else{
            return('NoArrayName')
        }
    }else{
        return('NoAffy')
    }
}

from_array_name_to_biomart_name = function(adf_file_dir,ebi_platform_name,biomart_array_name) {
    adf_file = paste0(adf_file_dir,ebi_platform_name,'.adf.txt')
    gse_array_name = catch_array_name(download_and_read_adf(adf_file,ebi_platform_name))
    if (!(grepl('No',gse_array_name))) {
        {
            if (sum(grepl(gse_array_name,biomart_array_name))) {
                biomart_name = paste(biomart_array_name[grep(gse_array_name,biomart_array_name)],collapse = ';')
                return(biomart_name)
            }

            else {
                return("Not Found")
            }
        }
    }
    else {
        return(gse_array_name)
    }
}


download_and_read_adf = function(adf_file_dir,ebi_platform_name) {
    adf_file = paste0(adf_file_dir,ebi_platform_name,'.adf.txt')
    if (!file.exists(adf_file)) {
        system(
            paste0(
                'wget https://www.ebi.ac.uk/arrayexpress/files/',as.character(ebi_platform_name),'/',
                as.character(ebi_platform_name),'.adf.txt ',
                '-O ',adf_file_dir,as.character(ebi_platform_name),'.adf.txt'
            )
        )
    }
    first_line_adf = catch_the_first_line_adf_file(adf_file)
    return(first_line_adf)
}




from_ebi_platform_name_to_biomart = function(ebi_platform_name,ebi_accession,gpl_df,biomart_array_name_list,adf_file_dir) {
    ebi_platform_name = as.character(ebi_platform_name)
    ebi_accession = as.character(ebi_accession)
    if (grepl('Illumina|Seq|need_test',ebi_platform_name)) {
        biomart_name = ebi_platform_name
    }else{
        if (grepl('GEOD',ebi_accession)) {
            gpl_name = from_ebi_GEOD_to_gpl(ebi_platform_name)
            biomart_name = from_gpl_id_to_gpl_title_then_biomart_name(gpl_name,gpl_df,biomart_array_name_list)
        }else{
            #adf_1st_line = read_sing_adf(adf_file_dir,ebi_platform_name)
            biomart_name = from_array_name_to_biomart_name(adf_file_dir,ebi_platform_name,biomart_array_name_list)
        }
    }
    return(biomart_name)
}


queryAE_info_pubmed = function(ebi_accession) {
    if (sum(grepl('Error',try(queryAE(keywords = ebi_accession))
    )) > 0) {
        pubmed = 'No'
    }
    else{
        temp_queryAE = queryAE(keywords = ebi_accession)
        file.remove(list.files(".", pattern = "\\.xml$"))

        pubmed = as.character(temp_queryAE$PubmedID)
    }
    return(pubmed)
}

remove_ebi_dataset_which_can_be_found_in_geo_search = function(ebi_summary_txt,geo_summary_xml) {
    row_logical = logical()
    for (i in 1:nrow(ebi_summary_txt)) {
        t_or_f = as.character(ebi_summary_txt$GSE[i]) %in% as.character(geo_summary_xml$GSE)
        row_logical = append(row_logical,t_or_f)
    }
    new_ebi_summary = ebi_summary_txt[!row_logical,]
    return(new_ebi_summary)
}

get_cols_contain_factor = function(df) {
    if (!as.logical(sum(grepl('LIBRARY_STRATEGY',colnames(df))))) {
        if (sum(grepl('^Array[.\\s]Data[.\\s]File',colnames(df)))) {
            col_name = (colnames(df)[grep('^Array[.\\s]Data[.\\s]File',colnames(df))])[1]
            if (length(df[,col_name]) == length(unique(as.character(df[,col_name])))) {
                if (sum(grepl(
                    '[Vv]alue|[Ff]actor|[tT]issue',colnames(df)
                )) > 0) {
                    df_contain_value = df[,c(col_name,colnames(df)[grep('[Ff]actor|[tT]issue',colnames(df))])]
                }else{
                    df_contain_value = as.data.frame(matrix(NA,0,ncol(df)))
                }
            }else {
                df_contain_value = 'repeated name'
            }
        }else{
            df_contain_value = as.data.frame(matrix(NA,0,ncol(df)))
        }
    }else {
        if (sum(grepl(
            'ENA_RUN|[Vv]alue|[Ff]actor|[tT]issue',colnames(df)
        )) > 0) {
            df_contain_value = df[,c("Source.Name",colnames(df)[grep('ENA_RUN|[Ff]actor|[tT]issue',colnames(df))])]
        }
        else {
            df_contain_value = as.data.frame(matrix(NA,0,ncol(df)))
        }
    }
    return (df_contain_value)
}

# Function
paste_df_in_one_row = function(df,sdrf_type) {
    if (is.data.frame(df)) {
        if (nrow(df) > 0) {
            new_row = character()
            for (i in 1:nrow(df)) {
                for (j in 1:(ncol(df) - 1)) {
                    new_row = paste0(new_row,df[i,j],'\t')

                }
                new_row = paste0(new_row,df[i,ncol(df)],'\n')
            }
        }
        else{
            if (sdrf_type == 'all') {
                new_row = 'no factor cols'
            }else{
                if (sdrf_type == 'related') {
                    new_row = '0 related'
                }else{
                    new_row = 'not sure'
                }
            }
        }
    }else{
        new_row = df
    }
    return (new_row)
}


get_seq_type_list = function(df) {
    type_seq = character()
    for (i in 1:nrow(df)) {
        if (df$sequence_type[i] == -1) {
            type_seq = append(type_seq,'M')
        }else{
            type_seq = append(type_seq,'R')
        }
    }
    return(type_seq)
}


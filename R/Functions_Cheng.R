# Functions used by Yiming

get_protein_coding<-function(df, spe){
    if (spe=="mouse"){
        ds_name<-"mmusculus_gene_ensembl"
    } else if (spe=="human"){
        ds_name<-"hsapiens_gene_ensembl"
    } else if (spe=="opossum"){
        ds_name<-"mdomestica_gene_ensembl"
    }
    ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=ds_name, host = "jul2015.archive.ensembl.org")
    gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), filters="ensembl_gene_id",
                      values=row.names(df), mart=ensembl_mart)

    gene_name$external_gene_name<-toupper(gene_name$external_gene_name)

    df$ensembl_gene_id<-row.names(df)
    #out_data<-merge(df, gene_name, by.x="ensembl_gene_id",by.y="ensembl_gene_id")
    out_data<-merge(df, gene_name, by="ensembl_gene_id")
    out_data<-subset(out_data, gene_biotype=="protein_coding" & length(external_gene_name)>1)
    out_data<-subset(out_data, select=-c(ensembl_gene_id, gene_biotype))

    out_data<-aggregate(.~external_gene_name, data=out_data, mean)
    row.names(out_data)<-out_data$external_gene_name
    out_data<-subset(out_data, select=-c(external_gene_name))

    return(out_data)
}

get_col_data<-function(sample_names){
    #spe_tissue=gsub("\\(r(\\d+)\\)$","\\1", sample_names)
    spe_tissue=gsub("\\(r(\\d+)\\)$","", sample_names)
    colData<-data.frame(do.call('rbind', strsplit(spe_tissue,"\\_")))
    row.names(colData)<-sample_names
    colnames(colData)<-c("Tissue", "Study")
    return(colData)
}

put_together_studies<-function(s_files, sum_cutoff){
    gx_data<-data.frame()

    for (t_s in s_files){
        s_data<-read.table(t_s, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
        s_data<-s_data[rowSums(s_data)>sum_cutoff,]

        if (dim(gx_data)[1]==0){
            gx_data<-s_data
        } else {
            gx_data<-merge(gx_data, s_data, by=0)
            rownames(gx_data)<-gx_data$Row.names
            gx_data<-subset(gx_data, select=-Row.names)
        }
    }
    return(gx_data)
}

MGF_MR_Marker<-function(project_dir){
    # project_dir<-"/data/home/cheng/Adipocyte/"
    # source(paste0(project_dir, "rCode/Merge_Data_Functions.R"))
    marker_dir<-paste0(project_dir, "MarkerPrediction")
    if (!dir.exists(marker_dir)){
        dir.create(marker_dir, recursive=TRUE)
    }

    # MicroArray data
    platforms<-"MicroArray"
    microarray_brw_file<-paste0(project_dir, platforms, "/BrW_Markers/rbatch.data_BrW.txt")
    microarray_data<-read.table(microarray_brw_file, sep="\t", check.names=FALSE, row.names=1)
    tissue_type<-str_extract(colnames(microarray_data), "BAT|WAT")
    colnames(microarray_data)<-tissue_type

    microarray.res.list <- getMarkerGenes(microarray_data, samples2compare="all", annotate=FALSE, score.cutoff=1)
    names(microarray.res.list)
    microarray.BAT<-microarray.res.list[["BAT_markers"]]
    microarray.BAT<-gsub(" ", "", microarray.BAT)
    microarray.WAT<-microarray.res.list[["WAT_markers"]]
    microarray.WAT<-gsub(" ", "", microarray.WAT)

    m.BAT<-str_split(microarray.BAT, ":", simplify = TRUE)
    df.m.BAT<-as.data.frame(m.BAT, stringsAsFactors=FALSE)
    rownames(df.m.BAT)<-df.m.BAT$V1
    df.m.BAT[,1]<-NULL
    colnames(df.m.BAT)<-"Score.Cutoff"

    m.WAT<-str_split(microarray.WAT, ":", simplify = TRUE)
    df.m.WAT<-as.data.frame(m.WAT, stringsAsFactors=FALSE)
    rownames(df.m.WAT)<-df.m.WAT$V1
    df.m.WAT[,1]<-NULL
    colnames(df.m.WAT)<-"Score.Cutoff"

    write.table(df.m.WAT, paste0(marker_dir, "/MicroArray.WAT.txt"), sep="\t", quote=FALSE)
    write.table(df.m.BAT, paste0(marker_dir, "/MicroArray.BAT.txt"), sep="\t", quote=FALSE)

    # RNASeq data
    platforms<-"RNASeq"
    rnaseq_brw_file<-paste0(project_dir, platforms, "/BrW_Markers/rbatch.data_BrW.txt")
    rnaseq_data<-read.table(rnaseq_brw_file, sep="\t", check.names=FALSE, row.names=1)
    tissue_type<-str_extract(colnames(rnaseq_data), "BAT|WAT")
    colnames(rnaseq_data)<-tissue_type

    rnaseq.res.list <- getMarkerGenes(rnaseq_data, samples2compare="all", annotate=FALSE, score.cutoff=1)
    names(rnaseq.res.list)
    rnaseq.BAT<-rnaseq.res.list[["BAT_markers"]]
    rnaseq.BAT<-gsub(" ","", rnaseq.BAT)
    rnaseq.WAT<-rnaseq.res.list[["WAT_markers"]]
    rnaseq.WAT<-gsub(" ", "", rnaseq.WAT)

    r.BAT<-str_split(rnaseq.BAT, ":", simplify = TRUE)
    df.r.BAT<-as.data.frame(r.BAT, stringsAsFactors=FALSE)
    rownames(df.r.BAT)<-df.r.BAT$V1
    df.r.BAT[,1]<-NULL
    colnames(df.r.BAT)<-"Score.Cutoff"

    r.WAT<-str_split(rnaseq.WAT, ":", simplify = TRUE)
    df.r.WAT<-as.data.frame(r.WAT, stringsAsFactors=FALSE)
    rownames(df.r.WAT)<-df.r.WAT$V1
    df.r.WAT[,1]<-NULL
    colnames(df.r.WAT)<-"Score.Cutoff"

    write.table(df.r.WAT, paste0(marker_dir, "/RNASeq.WAT.txt"), sep="\t", quote=FALSE)
    write.table(df.r.BAT, paste0(marker_dir, "/RNASeq.BAT.txt"), sep="\t", quote=FALSE)

    # Merge the markers from MicroArray and RNASeq

    df.BAT<-merge(df.r.BAT, df.m.BAT, by=0)
    rownames(df.BAT)<-df.BAT$Row.names
    df.BAT[,1]<-NULL
    colnames(df.BAT)<-c("Cutoff.RNASeq", "Cutoff.MicroArray")
    df.BAT$Cutoff.RNASeq<-as.numeric(df.BAT$Cutoff.RNASeq)
    df.BAT$Cutoff.MicroArray<-as.numeric(df.BAT$Cutoff.MicroArray)

    df.WAT<-merge(df.r.WAT, df.m.WAT, by=0)
    rownames(df.WAT)<-df.WAT$Row.names
    df.WAT[,1]<-NULL
    colnames(df.WAT)<-c("Cutoff.RNASeq", "Cutoff.MicroArray")
    df.WAT$Cutoff.RNASeq<-as.numeric(df.WAT$Cutoff.RNASeq)
    df.WAT$Cutoff.MicroArray<-as.numeric(df.WAT$Cutoff.MicroArray)

    df.BAT.gene<-get_protein_coding(df.BAT, "mouse")
    df.WAT.gene<-get_protein_coding(df.WAT, "mouse")
    df.BAT.gene<-df.BAT.gene[order(df.BAT.gene$Cutoff.RNASeq), ]
    df.WAT.gene<-df.WAT.gene[order(df.WAT.gene$Cutoff.RNASeq), ]

    write.table(df.WAT, paste0(marker_dir, "/marker.WAT.txt"), sep="\t", quote=FALSE)
    write.table(df.BAT, paste0(marker_dir, "/marker.BAT.txt"), sep="\t", quote=FALSE)
    write.table(df.WAT.gene, paste0(marker_dir, "/marker.gene.WAT.txt"), sep="\t", quote=FALSE)
    write.table(df.BAT.gene, paste0(marker_dir, "/marker.gene.BAT.txt"), sep="\t", quote=FALSE)

    write.table(rbind(df.BAT, df.WAT), paste0(marker_dir, "/marker.txt"), sep="\t", quote=FALSE)
    write.table(rbind(df.BAT.gene, df.WAT.gene), paste0(marker_dir, "/marker.gene.txt"), sep="\t", quote=FALSE)
}

Gen_ComData<-function(project_dir, marker_fd, spe_type){
    out_dir<-paste(project_dir, marker_fd, sep="/")
    if (!dir.exists(out_dir)){
        dir.create(out_dir, recursive=TRUE)
    }

    if (spe_type=="mouse"){
        # Retrieve only the datasets of BrW
        microarray_studies<-dir(paste0(project_dir, "MicroArray"), pattern="^[m]\\d+\\w_BrW", full.names=FALSE)
        # Remove m09M as pvWAT_M9 is clustered together with BAT, removed for marker prediction.
        microarray_studies<-microarray_studies[!grepl("m09M", microarray_studies)]

        study_files<-paste0(project_dir, "MicroArray","/", microarray_studies, "/", "ND_", microarray_studies, ".txt")
        microarray_all_data<-put_together_studies(study_files, 0)

        rnaseq_studies<-dir(paste0(project_dir, "RNASeq"), pattern="^[m]\\d+\\w_BrW", full.names=FALSE)
        study_files<-paste0(project_dir, "RNASeq","/", rnaseq_studies, "/tg/rlog.txt")
        rnaseq_all_data<-put_together_studies(study_files, -Inf)
        all_data<-merge(microarray_all_data, rnaseq_all_data, by=0)
        rownames(all_data)<-all_data$Row.names
        all_data<-all_data[,-1]
        batch_VM<-gsub(".*_([MR]\\d+)\\(r\\d+\\)", "\\1", colnames(all_data))
    } else if (spe_type=="human"){
        human_studies<-c("h01M", "h02M")
        study_files<-paste0(project_dir, "MicroArray","/", human_studies, "/", "ND_", human_studies, ".txt")
        all_data<-put_together_studies(study_files, 0)
        batch_VM<-strsplit2(colnames(all_data), split="_")[,2]
    }

    tissue_type<-str_extract(colnames(all_data), "BAT|WAT")
    tissue_type<-as.factor(tissue_type)

    if (length(unique(batch_VM))>1){
        gx_data_batch<-ComBat(dat = all_data, batch = batch_VM)
    } else {
        gx_data_batch<-all_data
    }
    write.table(all_data, paste0(out_dir, "/data_BrW.txt"), sep="\t", quote=FALSE)
    write.table(gx_data_batch, paste0(out_dir, "/rbatch.data_BrW.txt"), sep="\t", quote=FALSE)

    col_Data<-str_extract(colnames(gx_data_batch), "BAT|WAT")
    col_Data<-as.data.frame(col_Data, row.names=colnames(gx_data_batch))
    col_Data[,2]<-str_extract(colnames(gx_data_batch), "[MR]\\d+")
    colnames(col_Data)<-c("Tissue", "Study")

    if (length(unique(col_Data$Study))>1){
        design<-model.matrix(~0+Tissue+Study, data=col_Data)
    } else {
        design<-model.matrix(~0+Tissue, data=col_Data)
    }
    contrast_matrix = makeContrasts(TissueBAT-TissueWAT, levels=design)
    fit <- lmFit(all_data, design)
    fit_cont <- contrasts.fit(fit, contrast_matrix)
    fit_cont_eBayes <- eBayes(fit_cont)
    num_genes = dim(fit_cont_eBayes)[1]

    gene_list = topTable(fit_cont_eBayes, number=num_genes, sort.by="logFC")
    write.table(gene_list, file=paste0(out_dir, "/df_exp.txt"), sep="\t", quote=FALSE)

    gene_list = topTable(fit_cont_eBayes, number=num_genes, sort.by="logFC", lfc=1.5, p.value=1e-2, confint=TRUE)
    sig_gene_list<-get_protein_coding(gene_list, spe_type)

    write.table(sig_gene_list, file=paste0(out_dir, "/sig_gene.txt"), sep="\t", quote=FALSE)
}

MGF_MR_Marker_Com<-function(project_dir, marker_fd){
    # Use the combined data
    # Fold change cutoff
    fc_cutoff<-1.5

    # Load the data...
    # platform<-"BrW_Markers"
    brw_file<-paste0(project_dir, marker_fd, "/rbatch.data_BrW.txt")
    brw_data_all<-read.table(brw_file, sep="\t", check.names=FALSE, row.names=1)
    brw_data<-brw_data_all
    tissue_type<-str_extract(colnames(brw_data), "BAT|WAT")
    colnames(brw_data)<-tissue_type

    # Predict the markers genes...
    res.list <- getMarkerGenes(brw_data, samples2compare="all", annotate=FALSE, score.cutoff=1)
    names(res.list)
    all.BAT<-res.list[["BAT_markers"]]
    all.BAT<-gsub(" ", "", all.BAT)
    all.WAT<-res.list[["WAT_markers"]]
    all.WAT<-gsub(" ", "", all.WAT)

    BAT<-str_split(all.BAT, ":", simplify = TRUE)
    df.BAT<-as.data.frame(BAT, stringsAsFactors=FALSE)
    rownames(df.BAT)<-df.BAT$V1
    df.BAT[,1]<-NULL
    colnames(df.BAT)<-"Score.Cutoff"

    WAT<-str_split(all.WAT, ":", simplify = TRUE)
    df.WAT<-as.data.frame(WAT, stringsAsFactors=FALSE)
    rownames(df.WAT)<-df.WAT$V1
    df.WAT[,1]<-NULL
    colnames(df.WAT)<-"Score.Cutoff"

    # Load the differential gene expression
    df_exp<-read.table(paste0(project_dir, marker_fd, "/df_exp.txt"), sep="\t")

    # Filter out the markers based on FC for BAT
    markers.BAT<-rownames(df.BAT)
    df_exp_markers.BAT<-df_exp[markers.BAT,]
    df_exp_markers.BAT<-subset(df_exp_markers.BAT, abs(logFC)>fc_cutoff)
    df_exp_markers.BAT<-merge(df_exp_markers.BAT, df.BAT, by=0)
    rownames(df_exp_markers.BAT)<-df_exp_markers.BAT$Row.names
    df_exp_markers.BAT<-df_exp_markers.BAT[,-1]

    ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
    gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), filters="ensembl_gene_id",
                      values=row.names(df_exp_markers.BAT), mart=ensembl_mart)
    rownames(gene_name)<-gene_name$ensembl_gene_id
    df_exp_markers.BAT$GeneName<-gene_name[rownames(df_exp_markers.BAT), "external_gene_name"]
    df_exp_markers.BAT<-df_exp_markers.BAT[order(df_exp_markers.BAT$logFC),]
    write.table(df_exp_markers.BAT, paste0(project_dir, marker_fd, "/marker.BAT.", fc_cutoff, ".txt"), sep="\t", quote=FALSE)

    # Filter out the markers based on FC for WAT
    markers.WAT<-rownames(df.WAT)
    df_exp_markers.WAT<-df_exp[markers.WAT,]
    df_exp_markers.WAT<-subset(df_exp_markers.WAT, abs(logFC)>fc_cutoff)
    df_exp_markers.WAT<-merge(df_exp_markers.WAT, df.WAT, by=0)
    rownames(df_exp_markers.WAT)<-df_exp_markers.WAT$Row.names
    df_exp_markers.WAT<-df_exp_markers.WAT[,-1]

    gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), filters="ensembl_gene_id",
                      values=row.names(df_exp_markers.WAT), mart=ensembl_mart)
    rownames(gene_name)<-gene_name$ensembl_gene_id
    df_exp_markers.WAT$GeneName<-gene_name[rownames(df_exp_markers.WAT), "external_gene_name"]
    df_exp_markers.WAT<-df_exp_markers.WAT[order(df_exp_markers.WAT$logFC),]
    write.table(df_exp_markers.WAT, paste0(project_dir, marker_fd, "/marker.WAT.", fc_cutoff, ".txt"), sep="\t", quote=FALSE)

    write.table(rbind(df_exp_markers.BAT,df_exp_markers.WAT) , paste0(project_dir, marker_fd, "/marker.", fc_cutoff, ".txt"), sep="\t", quote=FALSE)

    # Combine WAT and BAT
    markers_data<-brw_data_all[c(rownames(df_exp_markers.BAT),rownames(df_exp_markers.WAT)), ]
    gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), filters="ensembl_gene_id",
                      values=row.names(markers_data), mart=ensembl_mart)
    rownames(gene_name)<-gene_name$ensembl_gene_id
    rownames(markers_data)<-gene_name[rownames(markers_data), "external_gene_name"]

    # Save the figures.
    pdf(file = paste0(project_dir, marker_fd, "/markers.heatmap.pdf"), onefile = FALSE)
    heatmap(as.matrix(markers_data), margins=c(5,10), Rowv=NA, cexCol=0.3, cexRow=0.5)
    dev.off()
}

pheatmap_plot<-function(test_data, test_fig_out,  use_heatmap=0, train_data=NULL, train_fig_out=NULL){
    # Plot the testing data ...
    test_data<-as.data.frame(t(test_data))
    r_names<-rownames(test_data)

    ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
    gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id",
                          values=r_names, mart=ensembl_mart)
    rownames(gene_name)<-gene_name$ensembl_gene_id

    #test_data<-test_data[gene_name$ensembl_gene_id, ]
    row.names(test_data)<-gene_name[r_names, "external_gene_name"]

    col_font_size<-ceiling(300/dim(test_data)[2])
    if (col_font_size>10){
        col_font_size<-10
    } else if (col_font_size<3){
        col_font_size<-3
    }

    f_width=dim(test_data)[1]/3
    if (f_width<5){
        f_width<-5
    } else if (f_width>10){
        f_width<-10
    }
    f_height=dim(test_data)[2]/2
    if (f_height<10){
        f_height<-10
    } else if (f_height>20){
        f_height<-20
    }

    setEPS()
    postscript(test_fig_out)

    #pdf(file = test_fig_out, onefile = FALSE, width=f_width, height=f_height)
    if (use_heatmap==0){
        pheatmap(as.matrix(test_data), clustering_distance_rows="euclidean", fontsize_row=3, fontsize_col=col_font_size)
    } else {
        #colfunc <- colorRampPalette(c("blue", "red"))
        colfunc<-colorRampPalette(c("green", "orange", "red"))
        cluster.col<-hcluster(t(test_data), method="euclidean")
        dd<-as.dendrogram(cluster.col)

        t_name<-colnames(test_data)
        col_weight<-rep(1, length(t_name))
        names(col_weight)<-t_name
        col_weight[names(col_weight)=="BAT"]<-100
        col_weight[names(col_weight)=="WAT"]<-100
        dd.reorder<-reorder(dd, col_weight)

        h_labels<-gsub(".eps$", ".labels.txt", test_fig_out)
        write.table(labels(dd.reorder), file=h_labels, sep="\t", quote=FALSE)
        #heatmap(as.matrix(test_data), margins=c(10,5), Rowv=NA, col = colfunc(20), cexRow=0.1+1/log10(dim(test_data)[1]+100), cexCol=0.1+1/log10(dim(test_data)[2]+10))
        heatmap.2(as.matrix(test_data), Colv=dd.reorder, Rowv=NA, col=colfunc(20), cexRow=0.05+1/log10(dim(test_data)[1]+100), cexCol=0.1+1/log10(dim(test_data)[2]+10),
                  scale="row", density.info="none", margins=c(10,5), trace="none", srtCol=45, dendrogram="column")
    }
    dev.off()
    if (!is.null(train_data)){
        # Plot the training data ...
        train_data<-as.data.frame(t(train_data))

        ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
        gene_name = getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id",
                          values=row.names(train_data), mart=ensembl_mart)
        train_data<-train_data[gene_name$ensembl_gene_id, ]
        row.names(train_data)<-gene_name$external_gene_name

        col_font_size<-ceiling(300/dim(train_data)[2])
        if (col_font_size>10){
            col_font_size<-10
        } else if (col_font_size<3){
            col_font_size<-3
        }

        pdf(file = train_fig_out, onefile = FALSE)
        pheatmap(as.matrix(train_data), clustering_distance_rows="euclidean", fontsize_row=3, fontsize_col=col_font_size)
        dev.off()

    }
}

########################################
# This function is for Webserver
########################################
norm_raw_TAR_data<-function(tar_file, chip_platform, spe){
    # tar_file<-"/data/home/cheng/test/GSE8044_RAW.tar"
    # chip_platform<-"affy_mouse430_2"
    # spe<-"mouse"

    # getAE("E-MTAB-758", type="raw"): this is to download the raw CEL file.
    #Unpack the CEL files
    cel_dir<-dirname(tar_file)
    untar(tar_file, exdir=cel_dir)
    cels = list.files(cel_dir, pattern = "cel$", ignore.case=TRUE)
    if (length(cels)>0){ # CEL file
        data_type="CEL"
    } else { # CEL.gz
        cels = list.files(cel_dir, pattern = "cel.gz$", ignore.case=TRUE)
        if (length(cels)>0){
            sapply(paste(cel_dir, cels, sep="/"), gunzip)
        } else {
            return(NA)
        }
    }

    cels_new = list.files(cel_dir, pattern = "CEL$", ignore.case=TRUE, full.names=TRUE)
    if (spe=="mouse"){
        ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl", host = "dec2015.archive.ensembl.org")
    } else if (spe=="human"){
        ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "dec2015.archive.ensembl.org")
    }

    ##################################################
    # READ and Process the data
    ##################################################
    if (grepl("_st_", chip_platform) | grepl("_hta_", chip_platform) |grepl(".db", chip_platform)){
        raw.data=read.celfiles(verbose=TRUE, filenames=cels_new) #From bioconductor
        data.rma.norm=oligo::rma(raw.data, normalize=TRUE)
    } else {
        raw.data=ReadAffy(verbose=TRUE, filenames=cels_new) #From bioconductor
        data.rma.norm=affy::rma(raw.data, normalize=TRUE)
    }

    #Get the important stuff out of the data - the expression estimates for each array
    rma_exp = exprs(data.rma.norm)

    #Format values to 5 decimal places
    rma_exp = round(rma_exp, digits=5)
    rma.avg = avereps(rma_exp) # By default, average based on row names, that is probes.
    rma.avg<-as.data.frame(rma.avg)
    rma.avg$ProbeNames<-rownames(rma.avg)

    probes=row.names(rma.avg)
    probes<-gsub("\\.\\d$", "", probes)
    row.names(rma.avg)<-probes
    if (grepl(".db", chip_platform)) {
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
    return(ID_data_x)
}

norm_raw_Seqcount<-function(count_file){
    gx_data<-read.table(count_file, header=TRUE, check.names=FALSE, row.names=1)
    gx_data<-gx_data[rowSums(gx_data)>1, ]
    sample_names<-colnames(gx_data)
    sample_names<-gsub("\\(r(\\d+)\\)$","", sample_names)

    col_Data<-as.data.frame(sample_names)
    row.names(col_Data)<-colnames(gx_data)
    colnames(col_Data)<-c("Samples")

    dds<-DESeqDataSetFromMatrix(countData=gx_data, colData=col_Data, design=~1)
    rld <- rlog(dds, blind=TRUE)
    rlog_data<-assay(rld)
    return(rlog_data)
}

#####################################################
# Functions for comparing different algorithms
#####################################################
cal_sen<-function(file_name, sheet_name, xlsx_file_name){
    res_data<-read.table(file_name, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    tissue_type<-row.names(res_data)
    tissue_type<-gsub("^\\w+\\.","",tissue_type)
    tissue_type<-gsub("\\.\\d","",tissue_type)
    res_data$tissue_type<-tissue_type

    sensi<-c()
    for (ii in 1:ncol(res_data)){
        sensi<-c(sensi, sum(res_data[,ii]==res_data$tissue_type)/dim(res_data)[1])
    }
    res_data["Accuracy",]<-as.character(sensi)

    write.xlsx(res_data, file=xlsx_file_name, sheetName=sheet_name,
               row.names=TRUE, append=TRUE, showNA=TRUE)

    res_data
}

cal_plot_sen<-function(base_dir, data_type, sig_gene_type, level_type){
    Sen_Excel_file<-paste0(base_dir, "ML/", level_type, ".Accuracy.xlsx")
    if (file.exists(Sen_Excel_file)){
        file.remove(Sen_Excel_file)
    }
    fig_dir<-paste0(base_dir, "ML/plots.", level_type, "/")

    sen_data<-data.frame()
    for (d_type in data_type){
        for (s_type in sig_gene_type){
            file_name<-paste0(base_dir, "ML/",level_type, ".", d_type, ".", s_type, ".txt")
            res_data<-cal_sen(file_name, paste(d_type, s_type, sep="."), Sen_Excel_file)
            res_data_sen<-res_data["Accuracy", ]
            rownames(res_data_sen)<-paste0(d_type,".",s_type)
            sen_data<-rbind(sen_data, res_data_sen)
        }
    }

    sen_data<-subset(sen_data, select=-c(tissue_type))
    ######################################
    # plot all the data in one Figure
    ######################################
    #par(mar=c(5,4,4,5), xpd=TRUE)
    par(mar=c(5,4,4,2), xpd=FALSE)
    matplot(t(sen_data), type="l", xlab="Machine Learning Algorithms", ylab="Accuracy", xaxt="n", ylim=c(0.3,1), xlim=c(1, ncol(sen_data)))
    matpoints(t(sen_data), type="p", pch=1:nrow(sen_data))
    axis(1, at = 1:ncol(sen_data), labels = colnames(sen_data), cex.axis =1)
    legend("topright", row.names(sen_data), pch=1:nrow(sen_data), cex=0.5)

    beige_acc<-31/37
    non_beige_RNASeq<-31/37
    non_beige_Microarray<-78/89
    abline(h=non_beige_RNASeq, lt=2) # RNAseq non-beige
    abline(h=non_beige_Microarray, lt=5) # Microarray non-beige
    text(4, non_beige_RNASeq, "non-beige RNASeq")
    text(4, non_beige_Microarray, "non-beige MicroArray")
    #png(file = "myplot.png", bg = "transparent")
    dev.copy(png,paste0(fig_dir, "accuracy.png"), width=960, height=960, res=100)
    dev.off()

    # rn_sen_data<-levels(as.factor(gsub("\\.\\w+$", "", row.names(sen_data))))
    ######################################
    # For each data type, plot the data
    ######################################
    for (temp_type in data_type){
        t_sen_data<-data.frame(t(sen_data))
        sub_sen_data<-data.frame(t(select(t_sen_data, starts_with(eval(temp_type)))))
        matplot(t(sub_sen_data), type="l", xlab="Machine Learning Algorithms", ylab="Accuracy", xaxt="n", ylim=c(0.3,1), xlim=c(1, ncol(sub_sen_data)))
        matpoints(t(sub_sen_data), type="p", pch=1:nrow(sub_sen_data))
        axis(1, at = 1:ncol(sub_sen_data), labels = colnames(sub_sen_data), cex.axis =1)
        legend("topright", row.names(sub_sen_data), pch=1:nrow(sub_sen_data), cex=1)
        dev.copy(png,paste0(fig_dir, "accuracy.", temp_type, ".png"), width=960, height=960, res=100)
        dev.off()

        sub_sen_data <- apply(sub_sen_data, 2, as.numeric)
        boxplot(sub_sen_data, main=temp_type)
        dev.copy(png,paste0(fig_dir, "algorithms.", temp_type, ".png"), width=960, height=960, res=100)
        dev.off()
    }

    m_sen_data <- apply(sen_data, 2, as.numeric)
    rownames(m_sen_data)<-rownames(sen_data)

    par(las = 1)
    boxplot(m_sen_data, main="All Data")
    dev.copy(png,paste0(fig_dir, "algorithms.png"), width=960, height=960, res=100)
    dev.off()

    par(mar=c(15,4,2,4)) # all axis labels horizontal
    x<-boxplot(t(m_sen_data), xlab="", xaxt = "n", main="All Data")
    text(1:nrow(sen_data), par("usr")[3], srt = 45, adj = 1, labels = rownames(m_sen_data), xpd = TRUE)
    dev.copy(png,paste0(fig_dir, "datasets.png"), width=960, height=960, res=100)
    dev.off()

    # plot the 3 datasets.
    df_data_type<-data.frame()
    for (temp_sig in data_type){
        temp_sig_data<-select(data.frame(t(sen_data)), starts_with(temp_sig))
        temp_sig_data<-as.data.frame(unlist(temp_sig_data, use.names=FALSE))
        names(temp_sig_data)<-temp_sig
        temp_sig_data<-data.frame(t(temp_sig_data))
        df_data_type<-rbind(df_data_type, temp_sig_data)
    }
    df_data_type<-apply(df_data_type, 1, as.numeric)
    x<-boxplot(df_data_type, xlab="", xaxt = "n", main="All Data")
    text(1:ncol(df_data_type), par("usr")[3], srt = 45, adj = 1, labels = colnames(df_data_type), xpd = TRUE)

    dev.copy(png,paste0(fig_dir, "datasets.only.png"), width=960, height=960, res=100)
    dev.off()

    # plot the markers effect.
    df_data_type<-data.frame()
    for (s_type in sig_gene_type){
        temp_sig_data<-select(data.frame(t(sen_data)), ends_with(s_type))
        temp_sig_data<-as.data.frame(unlist(temp_sig_data, use.names=FALSE))
        names(temp_sig_data)<-s_type
        temp_sig_data<-data.frame(t(temp_sig_data))
        df_data_type<-rbind(df_data_type, temp_sig_data)
    }
    df_data_type<-apply(df_data_type, 1, as.numeric)
    x<-boxplot(df_data_type, xlab="", xaxt = "n", main="All Data")
    text(1:ncol(df_data_type), par("usr")[3], srt = 45, adj = 1, labels = colnames(df_data_type), xpd = TRUE)

    dev.copy(png,paste0(fig_dir, "markers.only.png"), width=960, height=960, res=100)
    dev.off()

}

cal_plot_sen_BrBeW<-function(base_dir, level_type){
    Sen_Excel_file<-paste0(base_dir, "ML/", level_type, ".Accuracy.xlsx")
    if (!file.exists(Sen_Excel_file)){
        return()
    }
    fig_dir<-paste0(base_dir, "ML/plots.", level_type, "/")

    file <- system.file("tests", Sen_Excel_file, package = "xlsx")
    wb <- loadWorkbook(Sen_Excel_file)
    sheets <- getSheets(wb)

    BeBrW_All<-c("Brown", "Beige", "White")
    for (BeBrW in BeBrW_All){
        accuracy_BeBrW<-matrix(nrow = length(sheets),ncol=7)
        for (ii in 1:length(sheets)){
            sheet_ii<-read.xlsx2(Sen_Excel_file, ii, stringsAsFactors = FALSE)
            sheet_ii_name<-names(sheets[ii])

            sheet_ii<-subset(sheet_ii, tissue_type==BeBrW)
            sheet_ii<-sheet_ii[,-1]

            num_rows<-dim(sheet_ii)[1]
            num_cols<-dim(sheet_ii)[2]
            one_row<-c()
            for (jj in 1:(num_cols-1)){
                sheet_ii_jj<-sum(ifelse(sheet_ii[,jj]==sheet_ii$tissue_type,1,0))
                accu<-sheet_ii_jj/num_rows
                one_row<-c(one_row, accu)
            }
            accuracy_BeBrW[ii, ]<-one_row
        }
        rownames(accuracy_BeBrW)<-names(sheets)
        colnames(accuracy_BeBrW)<-c("svm","nB","lr","nnet","rF","boosting","rPart")

        boxplot(accuracy_BeBrW, main=BeBrW, ylim=c(0,1))
        dev.copy(png,paste0(fig_dir, BeBrW, ".png"), width=960, height=960, res=100)
        dev.off()
    }
}

get_study_str<-function(t_study){
    a<-strsplit(t_study, "")
    a<-a[[1]]
    if (a[2] == "0"){
        s_str<-paste0(a[4], a[3])
    } else {
        s_str<-paste0(a[4], a[2], a[3])
    }
    return(s_str)


}

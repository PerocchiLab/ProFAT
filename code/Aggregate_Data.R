# The function is used to remove the outliers and then aggregate the data
# Then, we cluster the data using three different metric distance.

library(amap)

# MicroArray
project_dir<-"/data/home/cheng/Adipocyte/"

#platforms<-c("MicroArray", "RNASeq")
platforms<-c("MicroArray")

for (p_form in platforms){
    sub_dir<-paste0(p_form, "/All_BeBrW")
    data_file<-"rbatch.All.ND.13.txt"
    if (p_form == "RNASeq"){
        sub_dir<-paste0(sub_dir, "/tg")
        data_file<-"rbatch.rlog.txt"
    }

    data_file<-paste0(project_dir, sub_dir,  "/", data_file)

    rbatch_nd_data<-read.table(data_file, sep="\t", check.names=FALSE, row.names=1)
    t_data<-as.data.frame(t(rbatch_nd_data))
    t_data$Samples<-gsub("\\(r\\d+\\)$", "", rownames(t_data))
    t_data_agg<-aggregate(.~Samples, data=t_data, mean)
    row.names(t_data_agg)<-t_data_agg[,1]
    t_data_agg<-t_data_agg[,-1]

    data_agg<-as.data.frame(t(t_data_agg))
    out_data_file<-gsub("txt$", "agg.txt", data_file)

    data_agg<-round(data_agg, digits=3)
    write.table(data_agg, out_data_file, sep="\t", quote=FALSE)

    #dist_metric<-c("euclidean", "spearman", "pearson")
    dist_metric<-c("euclidean")
    for (t_metric in dist_metric){
        fit <- hcluster(t(data_agg), method=t_metric) #

        pdf(file = paste0(gsub("txt$", "", out_data_file), "hcluster.", t_metric, ".pdf"))
        par(mar=c(4,2,2,6))
        plot(as.dendrogram(fit), xlab="Height", main=paste0(p_form, " (", t_metric, ")"), edgePar=list(col="blue", lwd=4), horiz=T)
        dev.off()
    }
}


# Plot the predictin accuracy.
# Load the data set information
library(ggplot2)
library(gplots)
library(biomaRt)
library(limma)
library(sva)

project_dir<-"/data/home/cheng/Adipocyte/"
source(paste0(project_dir, "rCode/Functions_Cheng.R"))

dataset_info<-paste0(project_dir, "DataInfo/All_Dataset.txt")
out_dir<-paste0(project_dir, "MLResults/")

all_datasets<-read.table(dataset_info, sep="\t", stringsAsFactors = FALSE)
colnames(all_datasets)<-all_datasets[1,]
all_datasets<-all_datasets[-1, ]

accu_file<-"all.pred.wt.null.nnet.txt"
studies<-unique(all_datasets$DataSetID)
all_accu_data<-data.frame()
all_raw_data<-data.frame()

for (t_study in studies){
    if (grepl("^[mh]\\d+M", t_study)){
        study_accu<-paste0(project_dir, "MicroArray/", t_study, "/", accu_file)
        study_data<-paste0(project_dir, "MicroArray/", t_study, "/", "ND_", t_study, ".txt")
        label_order<-paste0(project_dir, "MicroArray/", t_study, "/heatmap.labels.txt")
    } else if (grepl("^[mh]\\d+R", t_study)){
        study_accu<-paste0(project_dir, "RNASeq/", t_study, "/tg/", accu_file)
        study_data<-paste0(project_dir, "RNASeq/", t_study, "/tg/", "rlog.txt")
        label_order<-paste0(project_dir, "RNASeq/", t_study, "/tg/heatmap.labels.txt")
    }
    if (t_study == "m01R_BeBr"){
        stop_here=1
    }
    if (grepl("^m", t_study)){
        spe="mouse"
    } else if (grepl("^h", t_study)) {
        spe="human"
    }

    accu_data<-read.table(study_accu, check.names=FALSE)
    exp_data<-read.table(study_data, check.names=FALSE)
    exp_data<-get_protein_coding(exp_data, spe)

    sample_ids<-all_datasets[all_datasets$DataSetID==t_study, "StudyID"]
    sample_ids<-sample_ids[!grepl("pgWAT\\(T\\)_R1", sample_ids)] # Remove the outliesr.
    sample_ids<-sample_ids[!grepl("iWAT\\(T\\)_M11\\(r2\\)", sample_ids)] # Remove the outliesr.

    study_accu_data<-accu_data[sample_ids,]
    study_exp_data<-exp_data[, sample_ids]

    if (dim(all_accu_data)[1]==0){
        all_accu_data<-study_accu_data
        all_raw_data<-study_exp_data
    } else {
        all_accu_data<-rbind(all_accu_data, study_accu_data)
        all_raw_data<-merge(all_raw_data, study_exp_data, by=0)
        rownames(all_raw_data)<-all_raw_data$Row.names
        all_raw_data<-subset(all_raw_data, select=-Row.names)
    }

    num_samples<-dim(study_accu_data)[1]
    lab_size<-20/num_samples
    if (lab_size>1){
        lab_size<-1
    }

    #####################################
    # Barplot figure
    #####################################
    eps_file<-paste0(out_dir, "twobar.", t_study, ".", accu_file,".eps")
    #pdf(file = paste0(out_dir, "twobar.", t_study, ".", accu_file,".pdf"), onefile = FALSE)
    setEPS()
    postscript(eps_file)
    par(mar=c(10,5,2,2))
    # study_accu_data<-study_accu_data[order(study_accu_data$BAT),]

    l_order<-read.table(label_order, header=TRUE, stringsAsFactors = FALSE)
    l_order<-l_order[l_order$x != "BAT" & l_order$x !="WAT",]
    study_accu_data<-study_accu_data[l_order, ]

    mp<-barplot(as.matrix(t(study_accu_data)), args.legend = list(bty="n"), beside = FALSE, ylim=c(0,1), col=c("brown", "papayawhip"),
            las=2, cex.names=lab_size, axes = FALSE, axisnames = FALSE)
    text(mp, par("usr")[3], labels = rownames(study_accu_data), srt = 45,  adj=c(1.1, 1,1), xpd = TRUE, cex=lab_size)
    dev.off()

    #pdf(file = paste0(out_dir, "onebar.", t_study, ".", accu_file,".pdf"), onefile = FALSE)
    #par(mar=c(10,5,2,2), yaxs="i")
    #     barplot(as.matrix(t(study_accu_data)), args.legend = list(bty="n"), cex.names=lab_size, horiz = TRUE, las=2, col=c("brown", "white"),
#             space=-0.1, border=NA, xlim=c(0,1), xlab="Pro(BAT)")
#     box(lty = '1373', col = 'black')

    #dev.off()
    study_accu_data$TissueType<-gsub("\\(r\\d+\\)", "", rownames(study_accu_data))

    study_accu_data$pBAT<-study_accu_data$BAT
    study_accu_data<-study_accu_data[,c("pBAT", "TissueType")]

    num_tissue_type<-length(unique(study_accu_data$TissueType))
    # boxplot(pBAT ~ TissueType, data=study_accu_data)
    #####################################
    # boxplot
    #####################################
    f_size<-100/num_tissue_type
    if (f_size>10){
        f_size<-12
    } else if (f_size<5){
        f_size<-5
    }

    p<-ggplot(study_accu_data, aes(x=TissueType, y=pBAT))
    p<-p+ylim(0, 1)
    p<- p+theme(axis.text=element_text(size=f_size), axis.title=element_text(size=12,face="bold"), axis.text.x=element_text(angle=45, hjust=1)) # change text size
    p<- p+theme(axis.line = element_line(colour = "black"),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())
    p<-p+labs(x="Samples", y="Pro(BAT)")

    #p<-p+geom_boxplot()  # boxplot
    #p<-p+geom_point()
    #p<-p+geom_boxplot(size=1, aes(fill = TC_Type))  # boxplot
    p_jitter<-p+geom_jitter(position=position_jitter(width=.1, height=0.01), size=5, colour="red")  # add jitter

    pdf(file = paste0(out_dir, "jitter.", t_study, ".", accu_file,".pdf"), onefile = FALSE)
    print(p_jitter)
    dev.off()

    p_boxplot<-p+geom_boxplot()
    p_boxplot<-p_boxplot+geom_jitter(position=position_jitter(width=.1, height=0.01), size=5, colour="red")  # add jitter
    pdf(file = paste0(out_dir, "boxplot.", t_study, ".", accu_file,".pdf"), onefile = FALSE)
    print(p_boxplot)
    dev.off()

}
batch_VM<-strsplit2(colnames(all_raw_data), split="_")[,2]
batch_VM<-gsub("\\(r\\d+\\)", "", batch_VM)
gx_data_batch<-ComBat(dat = all_raw_data, batch = batch_VM)

ucp1_data<-as.data.frame(t(gx_data_batch["UCP1", ]))
all_accu_data_UCP1<-cbind(all_accu_data, ucp1_data)
sample_study<-strsplit2(rownames(all_accu_data_UCP1), "_")
all_accu_data_UCP1$Tissue<-gsub("\\.\\w+", "", sample_study[,1])
all_accu_data_UCP1$Tissue<-gsub("[a-z]+", "", all_accu_data_UCP1$Tissue)
all_accu_data_UCP1$Study<-gsub("\\(r\\d+\\)", "", sample_study[,2])
all_accu_data_UCP1$TColor<-"brown"
all_accu_data_UCP1[all_accu_data_UCP1$Tissue=="WAT", "TColor"]<-"papayawhip"
all_accu_data_UCP1[all_accu_data_UCP1$Tissue=="BAT(T)", "TColor"]<-"brown2"
all_accu_data_UCP1[all_accu_data_UCP1$Tissue=="WAT(T)", "TColor"]<-"orange"

# ggplot(all_accu_data_UCP1) + geom_point(aes(UCP1, BAT, col=Tissue), size=5)  +
#     scale_colour_manual(values = c("brown", "brown2", "papayawhip", "orange")) +
#     geom_smooth(method='lm',formula=y~x)
#
# ggplot(all_accu_data_UCP1,aes(UCP1, BAT))+stat_summary(fun.data=mean_cl_normal) +
#     geom_smooth(method='lm',formula=y~x)


all_accu_data_UCP1<-subset(all_accu_data_UCP1, Study!="R0" & Study!="R11")

g<-ggplot(all_accu_data_UCP1, aes(x=UCP1, y=BAT)) +
    geom_point(colour=all_accu_data_UCP1$TColor, size=5) +
    scale_colour_manual(values = c("brown", "brown2", "papayawhip", "orange")) +
    geom_smooth(method="loess") +
    labs(y="Pro(BAT)")
ggsave(paste0(out_dir, "UCP1_vs_ProBAT.eps"), device="eps", pch)
ggsave(paste0(out_dir, "UCP1_vs_ProBAT.pdf"))

setEPS()
postscript(paste0(out_dir, "UCP1_vs_ProBAT.eps"))
#pdf(file = paste0(out_dir, "UCP1_vs_ProBAT.pdf"), onefile = FALSE)
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
u <- par("usr") # The coordinates of the plot area
rect(u[1], u[3], u[2], u[4], border=NA, col='grey90')
par(new=TRUE)
plot(all_accu_data_UCP1$UCP1, all_accu_data_UCP1$BAT, col=all_accu_data_UCP1$TColor,pch=16, xlab="Ucp1 Relative Expression", ylab="Pro(BAT)")
legend("topleft", col=c("brown2", "brown", "papayawhip", "orange"), legend=c("BAT", "BAT(T)", "WAT", "WAT(T)"),pch=16)
dev.off()


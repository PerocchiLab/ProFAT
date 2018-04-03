# Pipeline

# Generate the DataInfo/Core_Dataset.txt. This txt file is a copy-paste from TableS1_Yiming.xlsx

# Generate the configuration file for MicroArray and RNASeq analysis.
# >gen_MR_Dataset.pl

# Run the MicroArray analysis. This function takes DataInfo/MicroArray_Dataset.txt and generates the normalized value for each study
# >source("MicroArray_Analysis.R)

# Run the RNASeq analysis. This function takes DataInfo/RNASeq_Dataset.txt and generates the fpkm, raw count, rlog, nc.
# >batch.RNASeq_Analysis.pl
# >RNASeq_Analysis.R
# >copy_qc_RNASeq.pl : copy the QC file to the folder.

# Merge all the data together. Put the different studies together and do the batch effect.
# >Merge_Data.R

# Aggregate the Data by taking the average of the biological replicates.
# >Aggregate_Data.R

# Do the quality check. Check the bioligcal replicates clustering.
# >pheatmap_distance_among_samples.R
# Remove samples iWAT(T)_M11(r2), pgWAT(T)_R1(r1) by modify the MicroArray_Dataset.txt and RNASeq_Dataset.txt. Rename them to the older
# one and remove the outlier samples. Rename the studies to OLD. Then repeat the whole process again (From Line 8)

# Marker Prediction. Predict the marker genes using MGFM/MGFR program.
# >Predict_Gene_Markers.R

# Superivsed machine learning to predict all the dataset.


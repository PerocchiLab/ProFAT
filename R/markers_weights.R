# jiang
# 2017.11.09
# plot relative importance for marker genes
library(clusterGeneration)
library(nnet)
library(ggplot2)
library(plyr)
library(hashmap)

adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
source(paste0(code_path,'gar_fun.r'))

###############################################
## set color according to ENSEMBLE ID
markerFolder = paste0(adipocyte_dir,'BrW_Markers/')
## get colors for markers
mouseMarker = read.delim(file = paste0(markerFolder,'marker.txt'),sep = '\t',
                         header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Ensembl = rownames(mouseMarker)
mouseMarkerBat = read.delim(file = paste0(markerFolder,'marker.BAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarkerWat = read.delim(file = paste0(markerFolder,'marker.WAT.1.5.txt'),sep = '\t',
                            header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mouseMarker$Type = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'BAT','WAT'))
mouseMarker$Color = sapply(mouseMarker$GeneName,function(x) ifelse(x %in% mouseMarkerBat$GeneName,'brown2','papayawhip'))
hashmap_type = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$Type)
hashmap_colors = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$Color)
hashmap_names = hashmap(keys = mouseMarker$Ensembl, values = mouseMarker$GeneName)

#######################################
## train neural network model
#######################################
## model from function nnet
mod1 = nnet(rand.vars, y, size=8, linout=T)


#################################################
#### get weights
###############################################
#use the gar.fun function 
weights = data.frame(gar.fun('y',mod1,bar.plot = FALSE))
relative_weight = weights[, c("inp.cont", "rel.imp")]
relative_weight = relative_weight[order(relative_weight$rel.imp, decreasing = FALSE),]
relative_weight$features = factor(row.names(relative_weight))
relative_weight$type = factor(hashmap_type[[relative_weight$features]])
relative_weight$colors = hashmap_colors[[relative_weight$features]]
relative_weight$gene = hashmap_names[[relative_weight$features]]


###############################################
## plot weights
###############################################

ggplot(relative_weight, aes(x = features, y = rel.imp, fill = type)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('brown','papayawhip')) +
  theme(legend.position="none") # remove legend

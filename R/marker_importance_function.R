
plot_marker_importance = function(training_data,save_file, markers_df){
  ## create some dictionary
  hashmap_type = hashmap(keys = markers_df$Ensembl, values = markers_df$Type)
  hashmap_colors = hashmap(keys = markers_df$Ensembl, values = markers_df$Color)
  hashmap_names = hashmap(keys = markers_df$Ensembl, values = markers_df$GeneName)

  
  ## train the model
  pred_method = 'nnet'
  pre_process<-NULL # NULL, scale, center. corr
  t_length<-3 # grid tune length.
  
  control =  trainControl(method="LOOCV", classProbs=TRUE)
  model1 =  train(TissueType ~., data=training_data, method=pred_method, preProcess=pre_process, trControl=control, tuneLength = t_length)
  
  # get variable importance
  weights = data.frame(gar.fun(training_data$TissueType, model1, bar.plot = FALSE))
  relative_weight = weights[, c("inp.cont", "rel.imp")]
  relative_weight = relative_weight[order(relative_weight$rel.imp, decreasing = FALSE),]
  relative_weight$features = factor(row.names(relative_weight)[order(relative_weight$rel.imp, decreasing = FALSE)], levels = row.names(relative_weight)[order(relative_weight$rel.imp, decreasing = FALSE)])
  relative_weight$type = factor(hashmap_type[[relative_weight$features]])
  relative_weight$colors = hashmap_colors[[relative_weight$features]]
  relative_weight$gene = hashmap_names[[relative_weight$features]]
  relative_weight$gene = factor(relative_weight$gene[order(relative_weight$rel.imp, decreasing = FALSE)], levels = relative_weight$gene[order(relative_weight$rel.imp, decreasing = FALSE)])
  
  pdf(file = save_file)
  ggplot(relative_weight, aes(x = gene, y = rel.imp, fill = type)) + 
    geom_bar(stat = 'identity') +
    #scale_fill_manual(values = c('brown','papayawhip')) +
    scale_fill_manual(values = c('cyan3','cyan3')) +
    theme(legend.position="none") + # remove legend 
    theme(axis.text.x = element_text(size = 10, angle = 45)) +
    coord_flip(ylim = c(-1,1)) +
    ylab('Relative Importance') +
    xlab('Marker Genes')
  dev.off()
}
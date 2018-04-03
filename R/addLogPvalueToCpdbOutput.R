library(xlsx)
adipocyte_dir = '/data/home/share/Projects/Adipocyte/'
#code_path = '/data/home/jiang/projects/scripts/adipocyte/20170126/'
code_path = '/data/home/share/Projects/Adipocyte/rCode/'
plots_dir = '/data/home/jiang/projects/adipocyte/test_plots/'
source(paste0(code_path,'functions_plots.R'))

markerFolder = paste0(adipocyte_dir,'BrW_Markers/')
## biomart get genes

#############################################################################
batGOcpdb = read.delim(file = paste0(markerFolder,'cpdb/batGo.tab'),
                       sep = '\t',header = TRUE,stringsAsFactors = FALSE)
batGOcpdb$`-log10(p-value)` = -log10(batGOcpdb$p.value)
write.xlsx(batGOcpdb[batGOcpdb$term_level==2,],
           file = paste0(markerFolder,'cpdb/cpdb_BatGo.xls'),
           sheetName = 'Level2',row.names = FALSE)
write.xlsx(batGOcpdb[batGOcpdb$term_level==3,],
           file = paste0(markerFolder,'cpdb/cpdb_BatGo.xls'),
           sheetName = 'Level3',row.names = FALSE,append = TRUE)
write.xlsx(batGOcpdb[batGOcpdb$term_level==4,],
           file = paste0(markerFolder,'cpdb/cpdb_BatGo.xls'),
           sheetName = 'Level4',row.names = FALSE,append = TRUE)

watGOcpdb = read.delim(file = paste0(markerFolder,'cpdb/watGo.tab'),
                       sep = '\t',header = TRUE,stringsAsFactors = FALSE)
watGOcpdb$`-log10(p-value)` = -log10(watGOcpdb$p.value)
write.xlsx(watGOcpdb[watGOcpdb$term_level==4,],
           file = paste0(markerFolder,'cpdb/cpdb_WatGo.xls'),
           sheetName = 'Level4',row.names = FALSE)
#######################################################################
dealwithCpdbPathway = function(cpdbPathwayFilePath){
  batPathwayCpdb = read.delim(file = cpdbPathwayFilePath,
                              sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  batPathwayCpdb = batPathwayCpdb[order(batPathwayCpdb$p.value),]
  batPathwayCpdb$`-log10(p-value)` = -log10(batPathwayCpdb$p.value)
  sources = unique(batPathwayCpdb$source)
  for(i in 1:length(sources)){
    if(i==1){
      write.xlsx(batPathwayCpdb[batPathwayCpdb$source==sources[i],],
                 file = gsub('tab$','xls',cpdbPathwayFilePath),
                 sheetName = sources[i],row.names = FALSE)
    }else{
      write.xlsx(batPathwayCpdb[batPathwayCpdb$source==sources[i],],
                 file = gsub('tab$','xls',cpdbPathwayFilePath),
                 sheetName = sources[i],row.names = FALSE,append = TRUE)
    }
  }
}

cpdbAllMarkersPathway = paste0(markerFolder,'cpdb/cpdbAllMarkersPathways.tab')
dealwithCpdbPathway(cpdbPathwayFilePath = cpdbAllMarkersPathway)
cpdbBatMarkersPathway = paste0(markerFolder,'cpdb/cpdbBatPathways.tab')
dealwithCpdbPathway(cpdbPathwayFilePath = cpdbBatMarkersPathway)

#######################################################################

dealwithCpdbGO = function(cpdbGOFilePath){
  allGOcpdb = read.delim(file = cpdbGOFilePath,
                         sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  allGOcpdb = allGOcpdb[order(allGOcpdb$p.value),]
  allGOcpdb$`-log10(p-value)` = -log10(allGOcpdb$p.value)
  tmpLevel = unique(allGOcpdb$term_level)
  tmpCat = unique(allGOcpdb$term_category)
  for(level_i in 1:length(tmpLevel)){
    for(cat_i in 1:length(tmpCat)){
      if(level_i==1&cat_i==1){
        if(sum((allGOcpdb$term_level==tmpLevel[level_i]) & 
               (allGOcpdb$term_category== tmpCat[cat_i]))>0){
          write.xlsx(allGOcpdb[((allGOcpdb$term_level==tmpLevel[level_i]) & 
                                  (allGOcpdb$term_category== tmpCat[cat_i])),],
                     file = gsub('tab$','xls',cpdbGOFilePath),
                     sheetName = paste0('Level',tmpLevel[level_i],tmpCat[cat_i]),
                     row.names = FALSE) 
        }
      }else{
        if(sum((allGOcpdb$term_level==tmpLevel[level_i]) &
               (allGOcpdb$term_category== tmpCat[cat_i]))>0){
          write.xlsx(allGOcpdb[((allGOcpdb$term_level==tmpLevel[level_i]) & 
                                  (allGOcpdb$term_category== tmpCat[cat_i])),],
                     file = gsub('tab$','xls',cpdbGOFilePath),
                     sheetName = paste0('Level',tmpLevel[level_i],tmpCat[cat_i]),
                     row.names = FALSE,append=TRUE)
        }
      }
    }
  }
}


cpdbAllMarkersGO= paste0(markerFolder,'cpdb/cpdbAllMarkersGo.tab')
dealwithCpdbGO(cpdbGOFilePath = cpdbAllMarkersGO)
cpdbBatMarkersGO = paste0(markerFolder,'cpdb/cpdbBatGo.tab')
dealwithCpdbGO(cpdbGOFilePath = cpdbBatMarkersGO)

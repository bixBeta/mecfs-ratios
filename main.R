dds.deseq.list = readRDS("~/Desktop/clusteRexplorer/pca-explorer/inf-ratios/INF-condition-day-dds.deseq.list.RDS")

# # deseq2 min rep Inf ------------------------------------------------------
# 
# preLim.colData = sv4.meta %>% rownames_to_column("temp") %>% 
#   select(orig.ident,5:10) %>%
#   unique()
# 
# preLim.colData$condition = paste0(preLim.colData$phenoGroup, "_", preLim.colData$day)
# 
# colData = preLim.colData[order(preLim.colData$orig.ident, colnames(pb.matrix[[1]])),]
# table(colData$orig.ident == colnames(pb.matrix[["20"]]))
# 
# pb.matrix = pb.matrix[-22]
# 
# dds.list = list()
# for (i in 1:length(pb.matrix)) {
#   
#   dds.list[[i]] <- DESeqDataSetFromMatrix(countData = pb.matrix[[i]],
#                                           colData = colData[colData$orig.ident %in% colnames(pb.matrix[[i]]),], design = ~ day)
#   names(dds.list)[[i]] <- paste0("dds__", names(pb.matrix)[[i]])
# }
# 
# dds.deseq.list = lapply(dds.list, function(x){
#   x = DESeq(x, minReplicatesForReplace = Inf)
# })
# 
# 
# saveRDS(dds.deseq.list, "INF-condition-day-dds.deseq.list.RDS")
# 



# filter norm counts ------------------------------------------------------

inf.norm.counts = lapply(dds.deseq.list, function(x){
  counts(x, normalized = T)
})



select.day = function(cluster){
  
  meta = as.data.frame(colData(cluster))
  
  meta.filter.day1 = meta %>% filter(day == "D1") %>% select(orig.ident)
  meta.filter.day2 = meta %>% filter(day == "D2") %>% select(orig.ident)
  
  return(list(d1 = unname(unlist(meta.filter.day1)), d2 = unname(unlist(meta.filter.day2))))
}

# Cluster 0 Analysis
ids = select.day(cluster = dds.deseq.list$dds__0)
day1.ids = ids$d1
day2.ids = ids$d2


# get.counts = function(counts){
#   
#   counts.day1 = as.data.frame(counts) %>% select(day1.ids)
#   counts.day2 = as.data.frame(counts) %>% select(day2.ids)
#   
#   return(list(day1.norm.counts = counts.day1, day2.norm.counts = counts.day2))
#   
# }

splitCounts = lapply(dds.deseq.list, function(x){
  
  metadata = as.data.frame(colData(x))
  metadata$concat = paste0(metadata$orig.ident,"_",metadata$ENID, "_",metadata$day)
  
  norm.counts = DESeq2::counts(x, normalized = T)
  
  # print(colnames(norm.counts))
  
  colnames(norm.counts) <- metadata$concat
  #print(colnames(norm.counts))
  
  days = select.day(x)
  
  
  # filtered.counts.day1 = counts.day1[rowSums(counts.day1 > 100 ) > 10 , ]
  # filtered.counts.day2 = counts.day2[rowSums(counts.day2 > 100 ) > 10 , ]
  
  filtered.counts = norm.counts[rowSums(norm.counts > 100 ) > 15 , ]
  
  # To DO ----------- look at every gene and make a wise decision
  # filtered.counts = filtered.counts + 1
  
  counts.day1 = as.data.frame(filtered.counts) %>% rownames_to_column("gene") %>% select(matches("D1"), "gene") %>% column_to_rownames("gene")
  counts.day2 = as.data.frame(filtered.counts) %>% rownames_to_column("gene") %>% select(matches("D2"), "gene") %>% column_to_rownames("gene")
  
  enids = unique(metadata$ENID)
  ratios = list()
  
  for (i in 1:length(enids)) {
    
    ratios[[i]] <- as.data.frame(counts.day2[,grep(enids[i], colnames(counts.day2))]  / counts.day1[,grep(enids[i], colnames(counts.day1))])
    
    names(ratios)[[i]] <- enids[i]
    
    rownames(ratios[[i]]) <- rownames(filtered.counts)
    
    try(colnames(ratios[[i]]) <- enids[i])
    
  }
  
  # try log2 transforming the ratios -- post filtering(zi)
  
  
  return(list(metadata = metadata,
              norm.counts= norm.counts, 
              days = days, 
              day1.norm.counts.filtered = counts.day1,
              day2.norm.counts.filtered = counts.day2,
              # filtered.counts.day1 = filtered.counts.day1,
              # filtered.counts.day2 = filtered.counts.day2 
              filtered.counts.100.15 = filtered.counts,
              ratios = ratios
  ))
  
})

library(purrr)

test.list = list()

for (i in 1:length(splitCounts)) {
  
  test.list[[i]]= pluck(splitCounts, i, "ratios")
  names(test.list)[[i]] <- names(splitCounts)[[i]]
}


ratios.list = list()

for (i in 1:length(splitCounts)) {
  
  ratios.list[[i]]= pluck(splitCounts, i, "ratios")
  names(ratios.list)[[i]] <- names(splitCounts)[[i]]
}

ratios.cbind = lapply(ratios.list, function(x){
  as.matrix(do.call(cbind, x))
})



clean = function(x){
  x[sapply(x, simplify = 'matrix', is.infinite)] <- 0
  x[sapply(x, simplify = 'matrix', is.nan)] <- 0
  x[sapply(x, simplify = 'matrix', is.na)] <- 0
}

x = lapply(ratios.cbind, function(i) if(is.numeric(i)) ifelse(is.infinite(i), 0, i) else i)
y = lapply(x, function(i) if(is.numeric(i)) ifelse(is.nan(i), 0, i) else i)
z = lapply(y, function(i) if(is.numeric(i)) ifelse(is.na(i), 0, i) else i)

getPCA = function(zi){
  
  
  rv <- rowVars(zi)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(zi[select,]))
  
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  pVar.df <- as.data.frame(percentVar)
  pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
  
  pVar.df = pVar.df[ , order(names(pVar.df))]
  pVar.df$percentVar = pVar.df$percentVar * 100
  pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
  
  
  return(list(scree = pVar.df, PCA = pca))
  
}

pcas = lapply(z, function(x){
  getPCA(zi = x)
})

saveRDS(pcas, "ratios.pcas.RDS")

pcas.df = lapply(pcas, function(x){
  as.data.frame(x[["PCA"]][["x"]][,1:2])
})




ggplotly(ggplot(pcas.df$dds__0 %>% rownames_to_column("ENID") %>% left_join(., colData, by = "ENID"),
                aes(x=PC1, y=PC2, colour = ENID)) + geom_point())

plot_ly(data = pcas$dds__2$scree, x = ~ x, y = ~ percentVar, type = "scatter", mode = "markers",
        width = 400, height = 400) %>%
  layout(autosize = F, 
         xaxis = list(categoryorder = "array",title = "nth PC", categoryarray = ~x),
         yaxis = list(title = "Percent Var", ticksuffix = "%"),
         title = "elbow plot")

# c0.d1 = splitCounts$dds__0$day1.norm.counts.filtered
# c0.d2 = splitCounts$dds__0$day2.norm.counts.filtered
# 
# enids = unique(splitCounts$dds__0$metadata$ENID)
# 
# d2.v.d1.matrix = counts.day2[,grep(enids[i], colnames(counts.day2))]  / counts.day1[,grep(enids[i], colnames(counts.day1))]
# 
# grep(pattern = enids[1], x = colnames(c0.d1))

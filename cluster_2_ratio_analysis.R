
# extract dds cluster2 ----------------------------------------------------

dds.c2 = dds.deseq.list$dds__2
colData.c2 = colData(dds.c2) %>% data.frame()
colData.c2$enid.id = paste0(colData.c2$orig.ident, "_", colData.c2$ENID)
norm.counts.c2 = counts(dds.c2, normalized = TRUE)

table(colnames(norm.counts.c2) == colData.c2$orig.ident)

colnames(norm.counts.c2) <- colData.c2$enid.id

# filter normalized counts to reduce noise --------------------------------

p6 = plotDensity(dds.list = dds.deseq.list, cluster = 2, jby = "orig.ident", minCounts = 20, minSamples = 15)

filtered.norm.counts.c2 = as.data.frame(p6$filteredCounts)
colnames(filtered.norm.counts.c2) = colnames(norm.counts.c2)

splitDays = function(filteredCounts, meta){
  day1.ids = meta %>% filter(day == "D1") %>% select(orig.ident, day)
  day2.ids = meta %>% filter(day == "D2") %>% select(orig.ident, day)
  
  counts.d1 = filteredCounts %>% select(matches(unlist(unname(day1.ids$orig.ident))))
  counts.d2 = filteredCounts %>% select(matches(unlist(unname(day2.ids$orig.ident))))
  
  return(list(
    filtered.counts.d1  = counts.d1,
    filtered.counts.d2  = counts.d2, 
    day1.ids = day1.ids,
    day2.ids = day2.ids
  ))
}

cluster2.split.counts = splitDays(filteredCounts = filtered.norm.counts.c2, meta = colData.c2)


# get ratio matrices ------------------------------------------------------

getRatio = function(splitCounts, colData){
  
  DAY1 = splitCounts[["filtered.counts.d1"]] + 0.1
  DAY2 = splitCounts[["filtered.counts.d2"]] + 0.1
  
  enids = unique(colData$ENID)
  ratios = list()
  
  for (i in 1:length(enids)) {
    
    ratios[[i]] <- as.data.frame(DAY2[,grep(enids[i], colnames(DAY2))]  / DAY1[,grep(enids[i], colnames(DAY1))])
    
    names(ratios)[[i]] <- enids[i]
    
    rownames(ratios[[i]]) <- rownames(DAY1)

    colnames(ratios[[i]]) <- enids[i]
  }
  
  return(ratios)
}

c2.ratios = getRatio(splitCounts = cluster2.split.counts, colData = colData.c2)
c2.ratios.cbind = as.matrix(do.call(cbind, c2.ratios))


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

c2.ratios.cbind.log10 = log10(c2.ratios.cbind)
c2.ratios.cbind.log2 = log2(c2.ratios.cbind)
c2.pca = getPCA(zi = c2.ratios.cbind.log10)
c2.pca.log2 = getPCA(zi = c2.ratios.cbind.log2)
ggplotly(ggplot(c2.pca$PCA$x[,c(1,2)] %>% data.frame() %>% rownames_to_column("ENID") %>% left_join(., colData.c2, by = "ENID") %>% filter(ENID !="E588") %>% filter(ENID != "E261") ,
                aes(x=PC1, y=PC2, colour = sex)) + geom_point())


plot_ly(data = c2.pca$scree, x = ~ x, y = ~ percentVar, type = "scatter", mode = "markers",
        width = 400, height = 400) %>%
  layout(autosize = F, 
         xaxis = list(categoryorder = "array",title = "nth PC", categoryarray = ~x),
         yaxis = list(title = "Percent Var", ticksuffix = "%"),
         title = "elbow plot")


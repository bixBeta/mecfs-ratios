dds.deseq.list = readRDS("~/Desktop/clusteRexplorer/pca-explorer/inf-ratios/INF-condition-day-dds.deseq.list.RDS")

source("../packages.R")
extract.filter.normCounts = function(dds){
  
  # extract coldata
  meta.data = colData(dds) %>% data.frame()
  
  # add enid.id column for grepping later
  meta.data$enid.id = paste0(meta.data$orig.ident, "_", meta.data$ENID)
  
  # extract normalized counts
  norm.counts = counts(dds, normalized = T) %>% data.frame()
  
  # add row medians to norm.counts
  norm.counts$row.medians = apply(X = norm.counts, MARGIN = 1, FUN = median)
  
  # get the quantile profile of the median
  q.tile = quantile(norm.counts$row.medians)
  

  q75.filtered.norm.counts = norm.counts[norm.counts$row.medians > q.tile[4], ]
  
  
  
  # plot norm counts pre and post filter
  
  counts.df.no.filter = norm.counts %>% select(-row.medians)
  counts.df.stacked.no.filter = stack(counts.df.no.filter)
  gg.df.no.filter = left_join(counts.df.stacked.no.filter, meta.data, by = c("ind" = "orig.ident"))
  
 p.pre.filter <- ggplot(gg.df.no.filter, aes(x=.data$values+1)) +
   stat_density(aes(group=.data$ind, color=.data$day), position="identity", geom="line", show.legend=TRUE) +
   scale_x_continuous(trans = log10_trans(),
                      breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(~10^.x))) +
   labs(color="") +
   xlab(paste0("normCounts_", deparse(substitute(dds)))) +
   ylab("Density") +
   ggtitle("Density of counts distribution") +
   theme_gray() + facet_wrap("day")


 counts.df.post.filter = q75.filtered.norm.counts %>% select(-row.medians)
 counts.df.stacked.post.filter = stack(counts.df.post.filter)
 gg.df.post.filter = left_join(counts.df.stacked.post.filter, meta.data, by = c("ind" = "orig.ident"))

 p.post.filter <- ggplot(gg.df.post.filter, aes(x=.data$values+1)) +
   stat_density(aes(group=.data$ind, color=.data$day), position="identity", geom="line", show.legend=TRUE) +
   scale_x_continuous(trans = log10_trans(),
                      breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(~10^.x))) +
   labs(color="") +
   xlab(paste0("normCounts_", deparse(substitute(dds)), "__", "Row_median >", round(q.tile[4], 2))) +
   ylab("Density") +
   ggtitle("Density of counts distribution") +
   theme_gray() + facet_wrap("day")


  
#return(p.pre.filter)
  return(list(colData = meta.data,
         norm.counts.No.filter = norm.counts,
         Quantile = q.tile,
         q75.filtered.norm.counts  = q75.filtered.norm.counts,
         ggplots = 
           list(ggplot.pre.filter = p.pre.filter,
           ggplot.post.filter = p.post.filter)

            )
         )
    
}

# x = extract.filter.normCounts(dds = dds.deseq.list$dds__19)

xl = lapply(dds.deseq.list, function(x){
  extract.filter.normCounts(dds = x)
})

saveRDS(xl, "filtered_counts_input_for_ratio_calculation.RDS")

density.plots.list = list()
for (i in 1:length(xl)) {
  density.plots.list[[i]] <- pluck(xl, i, "ggplots", 1) +  pluck(xl, i, "ggplots", 2) + ggtitle(names(xl)[[i]])
  names(density.plots.list)[[i]] <- names(xl)[[i]]
}

for (i in 1:length(density.plots.list)) {
  png(paste0("densityPlots/", names(density.plots.list)[[i]],".png"), width = 1920, height = 600, res = 150)
  print(density.plots.list[[i]])
  dev.off()
}

  # cluster 0 median histogram
c0.medians = xl[["dds__0"]][["norm.counts.No.filter"]]$row.medians %>% data.frame()
colnames(c0.medians) <- "values"

test <- ggplot(c0.medians , aes(x=.data$values+1)) +
  stat_density( position="identity", geom="line", show.legend=TRUE) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(~10^.x)))



# filtered.counts.list = list()
# for (i in 1:length(xl)) {
#   
#   filtered.counts.list[[i]] <-  pluck(xl, i, "q75.filtered.norm.counts")
#   names(filtered.counts.list)[[i]] <- names(xl)[[i]]
# }

top10 = paste0("dds__", 0:10)
top10.list = list()
for (i in 1:length(top10)) {
  
  top10.list[[i]] <- xl[[top10[i]]]
  names(top10.list)[[i]] <- top10[i]
 }



# calculate Ratios --------------------------------------------------------
# post extract.filter.normCounts func

getRatios = function(cluster){
  
  metadata = cluster[["colData"]] %>% data.frame()
  metadata$enid.id =  paste0(metadata$orig.ident,"_",metadata$ENID, "_",metadata$day)
  
  filtered.counts = cluster[["q75.filtered.norm.counts"]] %>% data.frame() %>% select(-row.medians)
  
  # return(filtered.counts %>% head()) 
  colnames(filtered.counts) <- metadata$enid.id
  
  meta.filter.day1 = metadata %>% filter(day == "D1") %>% select(orig.ident)
  meta.filter.day2 = metadata %>% filter(day == "D2") %>% select(orig.ident)
  
  d1 = unname(unlist(meta.filter.day1))
  d2 = unname(unlist(meta.filter.day2))
  
  counts.day1 = filtered.counts %>% rownames_to_column("gene") %>% select(matches("D1"), "gene") %>% column_to_rownames("gene")
  counts.day2 = filtered.counts %>% rownames_to_column("gene") %>% select(matches("D2"), "gene") %>% column_to_rownames("gene")
  
  counts.day1 = counts.day1 + 0.1
  counts.day2 = counts.day2 + 0.1
  
  enids = unique(metadata$ENID)
  ratios = list()
  
  for (i in 1:length(enids)) {
    
    ratios[[i]] <- as.data.frame(counts.day2[,grep(enids[i], colnames(counts.day2))]  / counts.day1[,grep(enids[i], colnames(counts.day1))])
    
    names(ratios)[[i]] <- enids[i]
    
    rownames(ratios[[i]]) <- rownames(counts.day2)
    
    colnames(ratios[[i]]) <- enids[i]
    
  }
  
  ratios.matrix =  as.matrix(do.call(cbind, ratios))
  log10.ratios.matrix = log10(ratios.matrix)
  
  return(list(D1 = counts.day1,
              D2 = counts.day2, 
              raw.ratios.list = ratios,
              non.log.ratios = ratios.matrix,
              log.10.ratios = log10.ratios.matrix,
              colData = metadata))
}
top10.Ratios = lapply(X = top10.list, function(x){
  getRatios(cluster = x)
})


counts.forPCA = lapply(top10.Ratios, function(x){
  chuck(x, "log.10.ratios")
})


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

pcas = lapply(counts.forPCA, function(x){
  getPCA(zi = x)
})

pca.plotly <- function(cluster, rm=NULL, pcx, pcy, overlay = "ENID"){
  
  #rm = c("E171", "E408" )
  
  pcas.subset.cluster.matrix = pcas[[paste0("dds__",cluster)]][["PCA"]][["x"]]
  pcas.scree = pcas[[paste0("dds__",cluster)]][["scree"]]
  
  meta.data = top10.Ratios[[paste0("dds__",cluster)]][["colData"]] %>% select(ENID, phenoGroup, sex) %>% unique()
  
  df = pcas.subset.cluster.matrix %>% data.frame() %>% rownames_to_column("ENID") %>% 
    left_join(.,meta.data, by = "ENID") %>% filter(!ENID %in% rm)
  
  u.pcx = paste0("PC", pcx)
  u.pcy = paste0("PC", pcy)
  
  p = ggplot(df ,
                  aes_string(x=df[,u.pcx], y=df[,u.pcy], colour = overlay, label = "ENID")) + geom_point() +
                 xlab(paste0(u.pcx, " : " , pcas.scree %>% filter(x == u.pcx) %>% select(percentVar) , "% variance")) +
                 ylab(paste0(u.pcy, " : ", pcas.scree %>% filter(x == u.pcy) %>% select(percentVar), "% variance")) +
                 coord_fixed()
               
               
               
  #print(p)
  return(p)
}

test = pca.plotly(cluster = 2, rm=NULL, pcx = 1, pcy = 3, overlay = "phenoGroup")
test

savePCAs = function(i){
  pc = pca.plotly(cluster = i, rm=NULL, pcx = 1, pcy = 2, overlay = "phenoGroup")
  pc2 = pca.plotly(cluster = i, rm=NULL, pcx = 1, pcy = 3, overlay = "phenoGroup")
  pc3 = pca.plotly(cluster = i, rm=NULL, pcx = 1, pcy = 4, overlay = "phenoGroup")
  pc4 = pc + pc2 + pc3
  return(pc4)
}
savePCAs(i = 0)



plotDensity = function(dds.list, cluster, jby = "orig.ident", minCounts, minSamples){
  
  dds = dds.list[[paste0("dds__",cluster)]]
  
  metadata = as.data.frame(colData(dds))
  metadata$concat = paste0(metadata$orig.ident,"_",metadata$ENID, "_",metadata$day)
  
  norm.counts = DESeq2::counts(dds, normalized = T)
  
  if (missing(minCounts)) {
    filtered.counts = norm.counts
  } else {
    filtered.counts = norm.counts[rowSums(norm.counts > minCounts ) > minSamples , ]
  }
  
  
  counts.df = as.data.frame(filtered.counts)
  counts.df.stacked = stack(counts.df)
  gg.df = left_join(counts.df.stacked, metadata, by = c("ind" = jby))
  
  p <- ggplot(gg.df, aes(x=.data$values+1)) +
    stat_density(aes(group=.data$ind, color=.data$day), position="identity", geom="line", show.legend=TRUE) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(~10^.x))) +
    labs(color="") +
    xlab(paste0("normCounts.", deparse(substitute(minCounts)), ".", deparse(substitute(minSamples)))) +
    ylab("Density") +
    ggtitle("Density of counts distribution") +
    theme_gray() + facet_wrap("day")
  
  return(list(densityPlot = p, filteredCounts = filtered.counts, args = c(minCounts = deparse(substitute(minCounts)), 
                                                                          minSamples = deparse(substitute(minSamples)))))
  
}
c19.m5 = cluster19.extractedCounts[["norm.counts.No.filter"]] %>% filter(row.medians > 5)

library(psych)
headTail(wc.input)
xtabs( ~ ENID + day,
       data = test)


gene.list = wc.input$geneID %>% unique()
wc.input.101 = wc.input %>% filter(geneID == "NOC2L")

wilcox.test(normCount ~ day,
            data = wc.input.101,
            paired = TRUE,
            conf.int = TRUE,
            conf.level = 0.95)


wc.test <- function(gene){
  
  wc.input.gene = wc.input %>% filter(geneID == gene)
  res = wilcox.test(normCount ~ day,
              data = wc.input.gene,
              paired = TRUE,
              conf.int = TRUE,
              conf.level = 0.95)
  
  return(list(res = res, wc.input = wc.input.gene))
  
}

case.Day2.wilcox = list()
for (i in 1:length(gene.list)) {
  
  case.Day2.wilcox[[i]] <- wc.test(gene = gene.list[i])
  names(case.Day2.wilcox)[[i]] <- as.character(gene.list[i])
  
}

case.day2.res = unlist(case.Day2.wilcox)
case.day2.pvals = as.data.frame(case.day2.res[grep(pattern = ".res.p.value", names(case.day2.res))])

case.day2.pvals$pvals = as.numeric(case.day2.pvals$`case.day2.res[grep(pattern = ".res.p.value", names(case.day2.res))]`)
case.day2.pvals$gene = rownames(case.day2.pvals)



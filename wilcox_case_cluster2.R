
# cluster2 filtered matrix  ----------------------------------------------

cluster2.filtered.matrix = xl[["dds__2"]][["q75.filtered.norm.counts"]] %>% data.frame() %>% select(-row.medians)
cluster2.metadata = xl$dds__2$colData


cluster2.cases.day1 = cluster2.metadata %>%
                            filter(phenoGroup == "G2" & day == "D1") %>% select(ENID, day) %>% rownames_to_column("orig.ident")
      
cluster2.cases.day2 = cluster2.metadata %>%
                            filter(phenoGroup == "G2" & day == "D2") %>% select(ENID, day) %>% rownames_to_column("orig.ident")

cluster2.cases.day1 = cluster2.cases.day1 %>% arrange(ENID)
cluster2.cases.day2 = cluster2.cases.day2 %>% arrange(ENID)

table(cluster2.cases.day1$ENID == cluster2.cases.day2$ENID)


transposed.c2.day1 = cluster2.filtered.matrix %>% select(cluster2.cases.day1$orig.ident)
transposed.c2.day2 = cluster2.filtered.matrix %>% select(cluster2.cases.day2$orig.ident)

colnames(transposed.c2.day1) <- cluster2.cases.day1$ENID
colnames(transposed.c2.day2) <- cluster2.cases.day2$ENID


transposed.c2.day1 = t(transposed.c2.day1)
transposed.c2.day2 = t(transposed.c2.day2)


c2.day1.matrix = melt(transposed.c2.day1)
c2.day2.matrix = melt(transposed.c2.day2)

c2.day1.matrix$day = "D1"
c2.day2.matrix$day = "D2"

colnames(c2.day1.matrix) <- c("ENID", "geneID", "normCount", "day" )
colnames(c2.day2.matrix) <- c("ENID", "geneID", "normCount", "day" )

wc.input = rbind(c2.day1.matrix, c2.day2.matrix)
wc.input = wc.input %>% select(geneID, day, ENID, normCount)

saveRDS(wc.input, "wilcoxon-case-input.RDS")

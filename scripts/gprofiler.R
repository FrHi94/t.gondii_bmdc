#Select Data from DGEA

#Set thresholds for the data in question:
de.markers <- subset(DEG12h, avg_logFC > 2 & p_val_adj < 0.01) 

#Do that for all the clusters in the DC data (including LPS):

#1.subset the data according to cluster identity 



library(gprofiler2)

gostres.BP <- gost(query = de.markers$external_gene_name, organism = "mmusculus", sources = "GO:BP")
#subset it with only one source e.g.: GO:BP! or Kegg etc.! 

#Plot the data
ggplot(head(gostres.BP$result, 20), aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value))) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Reds"))


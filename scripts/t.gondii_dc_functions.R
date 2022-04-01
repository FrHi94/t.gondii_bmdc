
##average heatmap function 

#data refers to a sc expression matrix you want to perform the analysis on 
#meta refers to the data you want to group your data by (e.g. clusters, celltypes, etc.)
#cond refers to the conditions you would like to display in the heatmap (numbers 1:16) (default is all conditions)
#seed refers to the features you would like to include (charactervector of genes to display)
#ann refers to the annotations you would like to have displayed in the heatmap  (default is all 3 annotations)

hmat_mousetoxo <- function(data,meta,cond = 1:16 , seed, ann = 1:3) {
  
  gene_data <- data.frame(data , cluster = meta ,check.names = F)  
  average_data <- aggregate(.~cluster, gene_data, mean)
  cluster_name <- average_data[,1]
  average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
  rownames(average_data) <- cluster_name #you can also give your data individual cluster names 
  average_data <- t(average_data)
  phmat1 <- t(scale(t(average_data)))
  #create matrix for averaged expression-values between -1.5 and 1.5 (this can be adjusted)
  phmat1[phmat1> 1.5] <- 1.5
  phmat1[phmat1 < -1.5] <- -1.5
  #plot genes of interest 
  phmat1 <- phmat1[,colorder] #insert if statement here to order the columns as you would like 
  
  pheatmap(phmat1[seed,cond], fontsize_row = 8, cellwidth = 15, clustering_distance_rows = "correlation", color = colorRampPalette(c("#440154" ,"#21908C", "#FDE725"))(200),
           annotation_col = annotation_col, annotation_colors = ann_cols,
           name = "value", column_split = annotation_col$CellType, 
           angle_col = "45"
  )
  
  
}


## Interaction plots

f.KNN_mousetoxo = function(seed, K=3, Nsteps=2, cutoffValue=0.65){
  
  # We first filter down the Top Interactors tables to the chosen number of interactors.
  TopIntGenes2 = TopIntGenes[1:K,]
  TopIntVal2 = TopIntVal[1:K,]
  
  # We then remove interactors with a weight below 0.65.
  if(cutoffValue<0){cutoffValue = 0; print("Warning: CutoffValue forced to 0!")}
  x.remove = abs(TopIntVal2) < cutoffValue
  TopIntGenes3 = TopIntGenes2; TopIntGenes3[x.remove] = NA
  TopIntVal3 = TopIntVal2; TopIntVal3[x.remove] = NA
  
  # Now we pick out the nearest neighbors, Nsteps steps away.
  seed2 = seed
  if(Nsteps>1){for(ii in 1:(Nsteps-1)){
    seed2 = c(seed2,TopIntGenes3[,seed2])
    seed2 = unique(seed2[!is.na(seed2)])
  }}
  TopIntGenes4 = TopIntGenes3[,seed2]
  TopIntVal4 = TopIntVal3[,seed2]
  
  # These objects are now transformed to an edgelist:
  keep.x = !is.na(as.vector(TopIntGenes4))
  elist1 = rep(colnames(TopIntGenes4), each=nrow(TopIntGenes4))[keep.x]
  elist2 = as.vector(TopIntGenes4)[keep.x]
  weights = as.vector(TopIntVal4)[keep.x]
  
  # There may be different fancy packages to create interactive plots, but I will use a 
  # straightforward solution: igraph
  
  gr.x = graph.data.frame(cbind(elist1,elist2), directed = T)
  
  # You can plot this graph, but it doesn't look very good. There are many things we can do to improve it. 
  V(gr.x)$label.family = "sans"
  V(gr.x)$label.cex = c(0.7)
  V(gr.x)$label.font = c(2)
  V(gr.x)$label.color = c("#333333")
  #V(gr.x)$label.dist = c(-0.5)
  #V(gr.x)$label.degree = pi/6
  E(gr.x)$arrow.size = 0.5
  E(gr.x)$width = abs(weights)*2
  V(gr.x)$name = make.names(V(gr.x)$name)
  #V(gr.x)$color[V(gr.x)$name %in% make.names(PTG_12h_markers)] = "#66CCFF"
  #V(gr.x)$color[V(gr.x)$name %in% make.names(LDM_12h_markers)] = "yellow3"
  V(gr.x)$color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_mouse))] = "#66CCFF"
  V(gr.x)$color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_toxo))] = "yellow3"
  V(gr.x)$frame.color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_mouse))] = "#3399CC"
  V(gr.x)$frame.color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_toxo))] = "yellow4"
  #V(gr.x)$shape[V(gr.x)$name %in% make.names(colnames(TopIntGenes_mouse))] = "circle"
  #V(gr.x)$shape[V(gr.x)$name %in% make.names(colnames(TopIntGenes_toxo))] = "circle"
  
  # We set the vertex sizes according to (minimum) distance from seed genes
  dist.x = distances(gr.x)
  for(ii in 1:ncol(dist.x)){dist.x[ii,ii] = NA}
  dist.x = dist.x[rownames(dist.x)%in%seed,]
  dist.x2 = apply(FUN=min, X=dist.x, MARGIN=2,na.rm=T)
  V(gr.x)$size = 2+10*(2^-dist.x2[V(gr.x)$name])
  
  
  # Color for interactions with positive correlations is set to red. 
  # Negative, if such exist, are colored blue.
  E(gr.x)$color[weights>0] = "#CC3333"
  E(gr.x)$color[weights<0] = "blue"
  
  # We can also map information to colors and shape of vertices.
  # There are two obvious options: 
  # 1) Different shapes for mouse and toxo, node colors according to a selected cluster from the umap.
  # 2) Make each node a pie chart with one piece for each cluster. Then distinction between mouse
  #    and toxo will have to be made some other way.
  #
  # I will leave this choice to a future discussion. For now, here is a simple function to translate
  # numeric values into a color scale, like a typical heatmap.
  
  f.TranslateToColor = function(val){
    bins.x = seq(from=-max(abs(val),na.rm=T)-0.01, to=max(abs(val),na.rm=T)+0.01, length.out=101)
    cval = c()
    for(ii in 1:length(val)){cval[ii] = sum(bins.x<val[ii])}
    cvec = colorRampPalette(c("blue", "dodgerblue", "white", "orange", "red"))(100)[cval]
    if(is.matrix(val)){
      cvec = matrix(cvec,nrow=nrow(val))
      rownames(cvec) = rownames(val)
      colnames(cvec) = colnames(val)
    } else {
      names(cvec) = names(val)
    }
    return(cvec)
  }
  par(mar=c(0,0,0,0)+.1)
  plot(gr.x, asp = 0) 
  return(gr.x)
  
  
}


##Subgraph of interaction plots

#gr = graphobject
#goi = gene of interest

subgraph_mtoxo = function(gr, goi){
  mem.x = clusters(gr)$membership
  index.x = mem.x[V(gr)$name == goi]
  gr.x2 = induced_subgraph(gr, which(mem.x==index.x))
  par(mar=c(0,0,0,0)+.1)
  plot(gr.x2, asp = 0) 
  return(gr.x2)
}


## Look at expression of cdks/cyclins across clusters /cell_class

#cluster = vector with cell and clusterannotation/celltype annotation or annotation of interest
#gene = vector with expression values for a gene of interest in the data
#dataset = which dataset to use (dc, dc+mac,mac)
#cols = colors used for the clusters 

CyCDK_expr <- function(object, clust, gene, cols){
  #cyclin B
  gg <- data.frame(cluster = clust, val = object@assays$SCT@data[gene,])
  # x should be the continuous variable and y the categorical group variable
  # You can color the ridges using any variable you want
  gg <- ggplot(gg, aes(x = val, y = cluster, fill = cluster)) +
    geom_density_ridges(alpha = 0.8) +# Set transparency
    scale_fill_manual(values = cols) +
    theme_ridges() +
    labs(x = "Expression value", y = "Cluster identity", fill = "Cluster identity") + 
    ggtitle(gene)
  return(gg)
}

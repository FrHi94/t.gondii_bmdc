hmat_mousetoxo <- function(data,meta) {
  
  
  gene_data <- data.frame(data , cluster = meta ,check.names = F)  
  average_data <- aggregate(.~cluster, gene_data, mean)
  cluster_name <- average_data[,1]
  average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
  rownames(average_data) <- cluster_name #you can also give your data individual cluster names 
  average_data <- t(average_data)
  phmat1 <- t(scale(t(average_data)))
  #create matrix for averaged expression-values between -1.5 and 1.5 (this can be adjusted) (clipping the values)
  #phmat1[phmat1> 1.5] <- 1.5
  # phmat1[phmat1 < -1.5] <- -1.5
  return(phmat1)
  
}



phmat_simple <- function(hmat, seed, cond) {
  pheatmap(hmat[seed,cond], fontsize_row = 4, cellwidth = 70,
           clustering_distance_rows = "correlation",breaks = seq(-max(abs(hmat)), max(abs(hmat)), length.out = 1000), 
           color = colorRampPalette(c("#660156","#21908C", "white", "#FFF68F", "#FFD700"))(1001)
        
           
  )
}
  
phmat_original <- function(hmat, seed, cond) {
  
  annotation_col <- data.frame(CellType = factor(rep(c("GM_DC", "GM_Mac"), each = 8)))
  rownames(annotation_col) <- colorder
  colorder <- c("ctrl_GM_DC", "PTG_Lys3h_GM_DC" , "PTG_3h_GM_DC", "PTG_12h_GM_DC", "LDM_Lys3h_GM_DC", "LDM_3h_GM_DC" , "LDM_12h_GM_DC", "LPS_GM_DC", "ctrl_GM_Mac", "PTG_Lys3h_GM_Mac" , "PTG_3h_GM_Mac", "PTG_12h_GM_Mac", "LDM_Lys3h_GM_Mac", "LDM_3h_GM_Mac" , "LDM_12h_GM_Mac", "LPS_GM_Mac")
  #add celltypes to the data in a complexheatmap
  annotation_col$Strain <- paste(c("Control", "PTG", "PTG", "PTG", "LDM", "LDM", "LDM", "LPS", "Control", "PTG", "PTG", "PTG", "LDM", "LDM", "LDM", "LPS"))
  annotation_col$Time <- as.integer(paste(c(3,3,3,12,3, 3, 12, 3, 3, 3, 3, 12, 3, 3,12, 3)))
  ann_cols <- list(CellType = c("GM_DC" = "mediumorchid4", "GM_Mac" = "turquoise4"), Strain = c("PTG" = "hotpink3", "LDM" = "olivedrab3", "Control" = "dodgerblue2", "LPS" = "red4"), Time = c("gray20", "gray87"))
  pheatmap(hmat[seed,cond], fontsize_row = 8, cellwidth = 15, clustering_distance_rows = "correlation",breaks = seq(-max(abs(hmat)), max(abs(hmat)), length.out = 1000), 
           color = colorRampPalette(c("#660156","#21908C", "white", "#FFF68F", "#FFD700"))(1001),
           annotation_col = annotation_col, annotation_colors = ann_cols,
           name = "value", column_split = annotation_col$CellType, 
           angle_col = "45"
  )
}


###co-expression network Toxoplasma 

f.KNN_toxo = function(seed, K=3, Nsteps=2, cutoffValue=0.65){
  
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
  E(gr.x)$width = abs(weights)*4
  V(gr.x)$name = make.names(V(gr.x)$name)
  #V(gr.x)$color[V(gr.x)$name %in% make.names(PTG_12h_markers)] = "#66CCFF"
  #V(gr.x)$color[V(gr.x)$name %in% make.names(LDM_12h_markers)] = "yellow3"
  V(gr.x)$color[V(gr.x)$name %in% make.names(crisp_INF_vs_naive$ensembl_gene_id)] = "#f4d35e" 
  #V(gr.x)$color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_toxo)[colnames(TopIntGenes_toxo) != crisp_INF_vs_naive$ensembl_gene_id])] = "grey"
  V(gr.x)$color[V(gr.x)$name %in% make.names(colnames(TopIntGenes_mouse))] = "#cc66cc"
  #V(gr.x)$shape[V(gr.x)$name %in% make.names(colnames(TopIntGenes_mouse))] = "circle"
  #V(gr.x)$shape[V(gr.x)$name %in% make.names(colnames(TopIntGenes_toxo))] = "circle"
  
  # We set the vertex sizes according to (minimum) distance from seed genes
  dist.x = distances(gr.x)
  for(ii in 1:ncol(dist.x)){dist.x[ii,ii] = NA}
  dist.x = dist.x[rownames(dist.x)%in%seed,]
  dist.x2 = apply(FUN=min, X=dist.x, MARGIN=2,na.rm=T)
  V(gr.x)$size = 2+5*(2^-dist.x2[V(gr.x)$name])
  
  
  # Color for interactions with positive correlations is set to red. 
  # Negative, if such exist, are colored blue.
  E(gr.x)$color[weights>0] = "#b86664"
  E(gr.x)$color[weights<0] = "grey"
  
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
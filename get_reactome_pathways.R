get_reactome_pathways <- function(){
  
  # get pathways from reactome
  pws = pathways("hsapiens", "reactome")
  
  # get pathways size
  pw.sizes = sapply(pws, function(x) x@protEdges[,c("src","dest")] %>% unlist %>% unique %>% length )
  
  # set minimum size to 10 genes
  pws = pws[pw.sizes > 10]
  pw.sizes = pw.sizes[pw.sizes > 10]
  
  # convert identifiers to symbol 
  pws = lapply(pws,convertIdentifiers, to = "symbol")
  
  # load STRING data 
  load("data/adj_string")
  string_genes = colnames(adj_string)
  
  # load mrna data 
  load("data/pancan_mrna_imputed_corrected_combined_not_reduced_genes")
  mrna_genes = lapply(mrna, colnames) %>% 
    unlist %>% table %>% `==`(12) %>% which %>% names
  
  # genes in both string and tcga all pancan12  
  common_genes = intersect(mrna_genes, string_genes)
  
  # genes in pws in string and in mrna  
  pws_genes = lapply(pws, function(x) intersect(
    x@protEdges[,c("src","dest")] %>% unlist %>% unique,
    common_genes))
  
  # calculate pathways coverages 
  smdf = data.frame( N0 = pw.sizes, 
                     N = sapply(pws_genes, length) )
  smdf$coverage = smdf$N/ smdf$N0
  
  smdf = smdf[order(smdf$coverage),]
  
  # pathways densities
  A = as.matrix(1*(adj_string != 0))
  # mean degrees
  pw_md = sapply(pws_genes, function(x){
    if(length(x) == 0 ) return(0)
    sum(rowSums(A[x,x,drop=FALSE]))/(2*length(x)^2)
  })
  smdf$md <- pw_md
  
  pathways_to_keep = rownames(smdf)[smdf$coverage > 0.5 & smdf$md>0 & smdf$N > 50 ]
  
  pws = pws[rev(pathways_to_keep)]
  pws_genes = pws_genes[rev(pathways_to_keep)]
  for(i in names(pws)){
    attr(pws[[i]], "genes_final") = pws_genes[[i]]
    attr(pws[[i]], "loginfo") = smdf[i,]
  }
  
  # save data to file
  SI = sessionInfo()
  save(file = "data/graphite_reactome_pathways", SI, pws)
}
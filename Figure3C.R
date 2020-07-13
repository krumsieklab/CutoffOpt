Figure3C <- function(cut_vec, data, adja, nboot, data_sizes){
  
  boots <- seq(from = 0, to = nboot, length = nboot+1)
  numb_edges <- length(which(adja == 1))
  
  cutoffs <- array(NaN, dim=c(length(data_sizes),length(cut_vec),length(boots)))
  fis_p_cut <- array(NaN, dim=c(length(data_sizes),length(cut_vec),length(boots)))
  
  index_data_size <- 0
  
  for (i_l in data_sizes) {
    
    index_data_size <- index_data_size + 1
    index_boot <- 0

    for (b in boots) {
      
      index_boot <- index_boot + 1
      index_cutoff <- 0

      if (b > 0) {
        data_re <- data[sample(nrow(data), i_l, replace = FALSE), ]   
      } else{
        data_re <- data[1:i_l,]
      }
      
      for (c in cut_vec) {
        
        index_cutoff <- index_cutoff + 1
        
        mc_tmp <- sco(data_re=data_re, i=i_l, cutoff=cut_vec[index_cutoff], adja=adja)
        
        cutoffs[index_data_size,index_cutoff,index_boot] <- mc_tmp[, 1]
        fis_p_cut[index_data_size,index_cutoff,index_boot] <- mc_tmp[, 2]
      }
    }
  }
  
  # compute confidence intervals across bootstrapping at the same size
  fis_cut_mean <- apply(fis_p_cut, c(1,2), function(x){
    confin(x)
  }) %>% .[2,,,drop=T]
  colnames(fis_cut_mean) <- format(round(cut_vec, 2), nsmall = 2)
  rownames(fis_cut_mean) <- data_sizes
  
  my_palette <- colorRampPalette(c("white", "black"))

  xlab <- colnames(fis_cut_mean)
  xlab[!((xlab %>% as.numeric) %in% as.character(seq(0,1,by = 0.1)))] <- ""
  
  # plot
  p <- pheatmap::pheatmap(-fis_cut_mean,col = my_palette(299), 
                     labels_col = xlab,
                     cluster_rows = F,cluster_cols = F,
                     show_rownames = T, show_colnames = T, main = "-Log10(Fisher's test p-value)")
  
  # save plot to file
  pdf("Figure3C.pdf", height = 15, width = 12)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  print(p)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Correlation cutoff", y=-0.02, gp=gpar(fontsize=11))
  grid.text("Sample size", x=0.97, rot=-90, gp=gpar(fontsize=11))
  dev.off()
  
}
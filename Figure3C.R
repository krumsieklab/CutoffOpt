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
  
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  lwid = c(0.2,4)
  lhei = c(0.8,4,0.8)
  
  pdf("Figure3C.pdf", height = 15, width = 12)
  gplots::heatmap.2(-fis_cut_mean, trace = "none", density.info = 'none', dendrogram = "none",
                    Rowv = NULL, Colv = NULL,
                    # #adjust breaks maybe
                    # breaks=seq(min(fis_cut_mean),0,length.out = 300),
                    col = my_palette(299), symkey = FALSE, keysize = 0.9,
                    # font size of rows and columns
                    cexRow=1.5,
                    cexCol=1.5,
                    ylab = "Sample size", xlab = "Correlation cutoff", main = "Fisher's test p-value (-log10)",
                    lmat = lmat, lwid = lwid, lhei = lhei
  )
  dev.off()
  
}
Figure4B <- function(cut_vec, data, adja, nboot, nswap){
  
  adja_once <- adja[1:20, 1:20]
  
  boots <- seq(from = 1, to = nboot, length = nboot)
  innerboots <- seq(from = 0, to = nboot, length = nboot+1)
  
  fis_p_swaps <- array(NaN, dim = c(length(nswap), length(cut_vec), length(boots), length(innerboots)))
  
  counter <- 0
  for(i in nswap) {
    
    counter <- counter+1
    
    for (b in boots) {
      
      # modify known adjancency matrix
      adja_rnd <- matrix(0L, nrow = dim(adja)[1], ncol = dim(adja)[2]) 
      adja_mod <- sym_generate_srand(adja_once, i) %>% as.matrix
      
      adja_rnd[1:20, 1:20] <- adja_mod
      adja_rnd[21:40, 21:40] <- adja_mod
      adja_rnd[41:50, 41:50] <- adja_mod[1:10, 1:10]
      
      for (ib in innerboots) {
        if (ib > 0) {
          data_re <- data[sample(nrow(data), nrow(data), replace = TRUE), ]
        } else{
          data_re <- data
        }
        
        gn_kor2 <- ggm.estimate.pcor(as.matrix(data_re), method = "dynamic", verbose=F)
        
        count_cutoff <- 0
        
        for (cutoff in cut_vec) {
          
          count_cutoff <- count_cutoff + 1
          adja_data <- gn_kor2
          
          adja_data[abs(adja_data) >= cutoff] <- 1
          adja_data[abs(adja_data) < cutoff] <- 0
          
          adja_data[lower.tri(adja_data, diag = TRUE)] <- NA
          adja_rnd[lower.tri(adja_rnd, diag = TRUE)] <- NA
          
          contin <- contab(adja_rnd, adja_data)
          
          # Fisher's exact test
          fis <- fisher.test(contin)
          fis_p_swaps[counter,count_cutoff,b,ib+1] <- -log10(fis$p.value)
          
        }
        
      }
      
    }
  }
  
  # aggregate and compute confidence intervals
  tmp <- apply(fis_p_swaps, c(1,2,3), mean)
  conf_fis <- apply(tmp, c(1,2), confin)
  
  rb <- palette(rainbow(length(nswap)-1))
  
  # merge data to plot in one data.frame
  dt <- lapply(1:dim(conf_fis)[2], function(x){
    cbind.data.frame(cut_vec=cut_vec, fis=conf_fis[2,x,],EdgeSwaps=sprintf("%d",nswap[x]))
  }) %>% do.call(rbind,.) %>% as.data.frame 
  
  # plot
  p <- ggplot(dt, aes(x = cut_vec, y = fis)) + 
    geom_line(aes(color = EdgeSwaps)) + 
    scale_color_manual(values = c("black", rb)) +
    theme_bw() +
    xlab("Correlation cutoff") +
    ylab("-log10(Fisher's test p-value)") +
    ggtitle("Incorrect Biological Reference")
  
  # save plot to file
  pdf("Figure4B.pdf", width = 10, height = 7)
  print(p)
  dev.off()
  
  return(p)
  
}
Figure4A <- function(cut_vec, data, adja, nboot, percentages){
  
  boots <- seq(from = 1, to = nboot, length = nboot)
  innerboots <- seq(from = 0, to = nboot, length = nboot+1)
  
  fis_p_allperc <- array(NaN, dim=c(length(percentages),length(cut_vec),length(boots)))
  
  counter <- 0
  
  for(percent in percentages){
    
    counter <- counter+1
    
    fis_p <- NULL
    count_cutoff <- 0
    
    for(cutoff in cut_vec){
      count_cutoff <- count_cutoff+1
      
      fis_tmp <- NULL
      adja_rnd <- adja
      
      for(b in boots) {
        fis_ib <- NULL
        
        #modify known adjancency matrix
        adja_mod <- adja[1:20,1:20]
        adja_mod[adja_mod == 1 & !is.na(adja_mod)] <- rand_adja(percent)
        
        adja_rnd[1:20, 1:20] <- adja_mod
        adja_rnd[21:40, 21:40] <- adja_mod
        adja_rnd[41:50, 41:50] <- adja_mod[1:10, 1:10]
        
        for(ib in innerboots){
          if (ib > 0) {
            data_re <- data[sample(nrow(data), nrow(data), replace = TRUE), ]
          } else{
            data_re <- data
          }
          
          #include Columns Age and Gender into calculation, but cut them out afterwards
          gn_kor2 <- ggm.estimate.pcor(as.matrix(data_re), method = "dynamic", verbose=F)
          
          gn_kor2[abs(gn_kor2) >= cutoff] <- 1
          gn_kor2[abs(gn_kor2) < cutoff] <- 0
          
          gn_kor2[lower.tri(gn_kor2, diag = TRUE)] <- NA
          adja_rnd[lower.tri(adja_rnd, diag = TRUE)] <- NA
          
          contin <- contab(adja_rnd, gn_kor2)
          
          #Fisher's exact test
          fis <- fisher.test(contin)
          fis_ib <- cbind(fis_ib, -log10(fis$p.value))
          
        }
        
        fis_tmp <- cbind(fis_tmp, mean(fis_ib))
        
      }
      
      fis_p_allperc[counter,count_cutoff,] <- fis_tmp
    }
    
  }
  
  # aggregate and compute confidence intervals
  tmp <- apply(fis_p_allperc, c(1,2,3), mean)
  conf_fis <- apply(tmp, c(1,2), confin)
  
  rb <- palette(rainbow(length(percentages)-1))
  
  # merge data to plot in one data.frame
  dt <- lapply(1:dim(conf_fis)[2], function(x){
    cbind.data.frame(cut_vec=cut_vec, fis=conf_fis[2,x,],MissingEdges=sprintf("%.0f%%",percentages[x]*100))
  }) %>% do.call(rbind,.) %>% as.data.frame 
  
  # plot
  p <- ggplot(dt, aes(x = cut_vec, y = fis)) + 
    geom_line(aes(color = MissingEdges)) + 
    scale_color_manual(values = c("black", rb)) +
    theme_bw() +
    xlab("Correlation cutoff") +
    ylab("-log10(Fisher's test p-value)") +
    ggtitle("Incomplete Biological Reference")
  
  # save plot to file
  pdf("Figure4A.pdf", width = 10, height = 7)
  print(p)
  dev.off()
  
  return(p)
  
}

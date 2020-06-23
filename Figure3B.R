Figure3B <- function(cut_vec, data, adja, nboot){
  
  fis_p <- NULL
  count_cutoffs <- 0
  
  for (cutoff in cut_vec) {
    fis_tmp <- cbind(cutoff)
    count_cutoffs <- count_cutoffs + 1
    
    for (i in 0:nboot) {
      if (i > 0) {
        data_re <- data[sample(nrow(data), nrow(data), replace = TRUE), ]
      } else{
        data_re <- data
      }
      
      # compute partial correlation matrix
      adja_data <-ggm.estimate.pcor(as.matrix(data_re), method = "dynamic", verbose=F)
      
      # binarize correlation matrix
      adja_data[abs(adja_data) >= cutoff] <- 1
      adja_data[abs(adja_data) < cutoff] <- 0
      
      adja[lower.tri(adja, diag = TRUE)] <- NA
      adja_data[lower.tri(adja_data, diag = TRUE)] <- NA
      
      # compute contingency table between GGM and prior knowledge
      contin <- contab(adja, adja_data)
      # compute fisher's test p-value
      fis <- fisher.test(contin)
      # store results
      fis_tmp <- cbind(fis_tmp, log10(fis$p.value))
      
    }
    fis_p <-rbind(fis_p, fis_tmp)
  }
  
  # compute confidence intervals
  conf <- t(apply(fis_p[,3:dim(fis_p)[2]], 1, confin))
  # merge with cutoff vector
  x <- cbind(fis_p[,1], conf)
  # assign names
  colnames(x) <- c("cut_vec","lowerbound","median_p","upperbound")
  
  p <- ggplot(as.data.frame(x), aes(x=cut_vec, y=-median_p))+
    geom_errorbar(aes(ymin=-upperbound, ymax=-lowerbound), colour="grey")+
    geom_line(size=1)+
    xlab("Correlation cutoff")+
    ylab("Fisher's test p-value (-log10)")+
    theme_light(base_size = 15) +
    ggtitle("GeneNet Cutoff optimization")
  
  pdf("Figure3B.pdf", width = 10.5, height= 10)
  print(p)
  dev.off()
  
}
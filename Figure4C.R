Figure4C <- function(cut_vec, data, adja, adja_block, adja_1s, nboot){
  
  fis_p <- NULL
  fis_p_block <- NULL
  fis_p_1s <- NULL
  count_cutoffs <- 0
  
  for (cutoff in cut_vec) {
    fis_tmp <- cbind(cutoff)
    fis_tmp_block <- cbind(cutoff)
    fis_tmp_1s <- cbind(cutoff)
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
      
      adja_data[lower.tri(adja_data, diag = TRUE)] <- NA
      adja[lower.tri(adja, diag = TRUE)] <- NA
      adja_block[lower.tri(adja_block, diag = TRUE)] <- NA
      adja_1s[lower.tri(adja_1s, diag = TRUE)] <- NA
      
      # compute contingency table between GGM and prior knowledge
      contin <- contab(adja, adja_data)
      contin_block <- contab(adja_block, adja_data)
      contin_1s <- contab(adja_1s, adja_data)
      # store results
      fis_tmp <- cbind(fis_tmp, log10(fisher.test(contin)$p.value))
      fis_tmp_block <- cbind(fis_tmp_block, log10(fisher.test(contin_block)$p.value))
      fis_tmp_1s <- cbind(fis_tmp_1s, log10(fisher.test(contin_1s)$p.value))
      
    }
    fis_p <-rbind(fis_p, fis_tmp)
    fis_p_block <-rbind(fis_p_block, fis_tmp_block)
    fis_p_1s <-rbind(fis_p_1s, fis_tmp_1s)
  }
  
  # compute confidence intervals
  conf <- t(apply(fis_p[,3:dim(fis_p)[2]], 1, confin))
  conf_block <- t(apply(fis_p_block[,3:dim(fis_p_block)[2]], 1, confin))
  conf_1s <- t(apply(fis_p_1s[,3:dim(fis_p_1s)[2]], 1, confin))
  # merge with cutoff vector
  x <- rbind.data.frame(cbind.data.frame(cut_vec=fis_p[,1], conf, Adjacency=rep("Original",dim(conf)[1])), 
                        cbind.data.frame(cut_vec=fis_p[,1], conf_block, Adjacency=rep("Subclass",dim(conf)[1])), 
                        cbind.data.frame(cut_vec=fis_p[,1], conf_1s, Adjacency=rep("1-Sugar",dim(conf)[1])))
  # assign names
  colnames(x) <- c("cut_vec","lowerbound","median_p","upperbound","Adjacency")
  
  p <- ggplot(as.data.frame(x), aes(x=cut_vec, y=-median_p, color=Adjacency))+
    geom_errorbar(aes(ymin=-upperbound, ymax=-lowerbound), colour="grey")+
    geom_line(size=1)+
    xlab("Correlation cutoff")+
    ylab("Fisher's test p-value (-log10)")+
    scale_colour_manual(values = c("Original"= "black","1-Sugar"="dodgerblue3", "Subclass"="seagreen3")) +
    theme_light(base_size = 15) +
    ggtitle("GeneNet Cutoff optimization")
  
  pdf("Figure4C.pdf", width = 10.5, height= 10)
  print(p)
  dev.off()
  
}
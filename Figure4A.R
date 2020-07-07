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

tmp <- apply(fis_p_allperc, c(1,2,3), mean)
conf_fis <- apply(tmp, c(1,2), confin)

#save max and min for different axis
maxYleft = 60
minYleft = 0

rb <- palette(rainbow(9))

pdf("Figure4A.pdf", width = 10, height = 7)
par(mar=c(5,4,4,6)+0.1)
plot(cut_vec, conf_fis[2,1,], type='l', axes=TRUE, col="black", xlab = "", ylab = "",
     ylim = c(minYleft, maxYleft))
# axis(2, ylim=c(minYleft, maxYleft), col="black", las=1)
mtext("Fisher's test p-value (log10)", side=2, line=2.5)
abline(h = 0, lty = 2, col = "black")
for(i in 2:dim(conf_fis)[2]){
  lines(cut_vec, conf_fis[2,i,], type="l", col=rb[i-1])
  # par(new=TRUE) 
  # plot(cut_vec, conf_fis[2,i,], type='l', axes=FALSE, col=rb[i-1], xlab="", ylab="",
  #      ylim = c(minYleft, maxYleft))
}
axis(1,xlim=c(0,1))
mtext("Correlation cutoff",side=1,col="black",line=2.5) 
legend(
  "bottomright",
  title = "Missing edges of original adjacency matrix",
  c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%"),
  lty = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  col = c("black", rb)
) 
dev.off()

}

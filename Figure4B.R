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

tmp <- apply(fis_p_swaps, c(1,2,3), mean)
conf_fis <- apply(tmp, c(1,2), confin)

#save max and min for different axis
maxYleft <- 60
minYleft <- 0

rb <- palette(rainbow(dim(fis_p_swaps)[1]+1))

pdf("Figure4B.pdf", width = 8, height = 7)
plot(cut_vec, conf_fis[2,1,], type='l', axes=TRUE, col="black", xlab = "", ylab = "",
     ylim = c(minYleft, maxYleft))
# axis(2, ylim=c(minYleft, maxYleft), col="black", las=1)
mtext("Fisher's test p-value (log10)", side=2, line=2.5)
for(i in 2:dim(fis_p_swaps)[1]){
  lines(cut_vec, conf_fis[2,i,], type="l", col=rb[i])
}
abline(h = 0, lty = 2, col = "black")
axis(1,xlim=c(0,1))
mtext("Correlation cutoff",side=1,col="black",line=2.5) 
legend(
  "topright",
  title = "Swapped Edges",
  nswap %>% as.character(),
  lty = rep(1,length(nswap)),
  col=c("black",rb[1:length(nswap)])
) 
dev.off()

}
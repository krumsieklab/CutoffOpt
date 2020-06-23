# Downloads files from the web and verifies their checksum. Will use local copy in current directory, if it exists
load.web.file <- function(
  url, md5sum, outfile, zipfile = F) {
  # check if local file exists
  if (file.exists(outfile)) {
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) stop(sprintf("Local file %s has wrong checksum: %s", outfile, realsum))
    # do not delete wrong file, it was already here before
    
  } else {
    if(zipfile){
      # download file
      temp <- tempfile()
      download.file(url,temp)
      unzip(zipfile = temp, files = outfile, exdir = ".")
    } else {
      # download file
      download.file(url, outfile)
    }
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) { 
      # delete wrong file
      unlink(outfile)
      stop(sprintf("Remote file %s has wrong checksum: %s", url, realsum))
    }
  }
}

# build contingency table from two adjacency matrices
contab <- function(adja, data) {
  ct_a <- length(which(adja == 1 & data == 1)) # true positive, upper left
  ct_b <- length(which(adja == 0 & data == 1)) # false positive, upper right
  ct_c <- length(which(adja == 1 & data == 0)) # false negative, lower left
  ct_d <- length(which(adja == 0 & data == 0)) # true negative, lower right
  ct <- cbind(c(ct_a, ct_c), c(ct_b, ct_d))
  return(ct)
}

#calculate 95% confidence interval of the values (not CI of mean)
confin <- function(fisp) {
  tmean <- mean(fisp,na.rm = T)
  tupper <- quantile(fisp, probs = (0.025), na.rm = T)
  tlower <- quantile(fisp, probs = (0.975), na.rm = T)
  conf <- cbind(tupper, tmean, tlower)
}

# create randomized adjacency matrix
rand_adja <- function(perc) {
  # in the original adjacency (per subclass) there are 27 edges
  # randomly remove perc% of them
  x27 <- runif(27, 0, 1)
  x27[x27 > perc] <- 1
  x27[x27 <= perc] <- 0
  return(x27)
}

sco <- function(data_re,i,cutoff,adja){
  
  # compute partial correlation matrix
  adja_data <- ggm.estimate.pcor(as.matrix(data_re), method = "dynamic", verbose=F)
  
  # compute adjacency matrix according to cutoff
  adja_edge <- adja_data
  adja_edge[abs(adja_edge) >= cutoff] <- 1
  adja_edge[abs(adja_edge) < cutoff] <- 0
  
  adja[lower.tri(adja, diag = TRUE)] <- NA
  adja_edge[lower.tri(adja_edge, diag = TRUE)] <- NA
  # compute contingency table
  contin_edge <- contab(adja, adja_edge)
  # compute Fisher's exact test at given cutoff
  fis_edge <- fisher.test(contin_edge) %>% .$p.value %>% log10
  
  # compute adjacency matrix according to a statistical cutoff of 0.01 FDR
  gn_pvalues <- network.test.edges(adja_data, plot=F, verbose=F)
  gn_pvalues$pval <- p.adjust(gn_pvalues$pval,"fdr")
  gn_pvalues$pval[gn_pvalues$pval > 0.01] <- NaN
  fdr_0.01 <- abs(gn_pvalues$pcor[which.max(gn_pvalues$pval)])
  if(length(fdr_0.01)==0){
    fdr_0.01 <- NaN
  }
  
  cbind(fdr_0.01,fis_edge)
  
}

# generate adjacency with ntry swaps preserving node degree
sym_generate_srand <- function(s1,ntry){
  
  nrew <- 0
  srand <- s1
  #save indices of fields that are 1 and only take upper right triangle
  index_srand <- which(srand != 0, arr.ind=T)
  index_srand <- index_srand[index_srand[,2] > index_srand[,1],]
  i_srand <- index_srand[,1]
  j_srand <- index_srand[,2]
  Ne <- length(i_srand)
  
  if(ntry==0) {
    srand=s1
  } else {
    for (i in 1:ntry){
      e1 <- 1+floor(Ne*runif(1, 0, 1))
      e2 <- 1+floor(Ne*runif(1, 0, 1))
      v1 <- i_srand[e1]
      v2 <- j_srand[e1]
      v3 <- i_srand[e2]
      v4 <- j_srand[e2]
      if ((v1!=v3)&&(v1!=v4)&&(v2!=v4)&&(v2!=v3)){
        if (runif(1, 0, 1)>0.5){
          if ((srand[v1,v3]==0)&&(srand[v2,v4]==0)){
            srand[v1,v2] <- 0
            srand[v3,v4] <- 0
            srand[v2,v1] <- 0
            srand[v4,v3] <- 0
            srand[v1,v3] <- 1
            srand[v2,v4] <- 1
            srand[v3,v1] <- 1
            srand[v4,v2] <- 1
            nrew <- nrew+1
            i_srand[e1] <- v1
            j_srand[e1] <- v3
            i_srand[e2] <- v2
            j_srand[e2] <- v4
          }
        } else {
          v5 <- v3
          v3 <- v4
          v4 <- v5
          rm(v5)
          if ((srand[v1,v3]==0)&&(srand[v2,v4]==0)){
            srand[v1,v2] <- 0
            srand[v4,v3] <- 0
            srand[v2,v1] <- 0
            srand[v3,v4] <- 0
            srand[v1,v3] <- 1
            srand[v2,v4] <- 1
            srand[v3,v1] <- 1
            srand[v4,v2] <- 1
            nrew <- nrew+1
            i_srand[e1] <- v1
            j_srand[e1] <- v3
            i_srand[e2] <- v2
            j_srand[e2] <- v4
          }
        }
      }
    }
  }
  return(srand)
}


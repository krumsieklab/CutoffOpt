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

# Adaptation of GeneNet (version  1.2.15) package code to run silently
# source code taken directly from package
# (removed hardcoded, unconditional calls of cat function)
network.test.edges.silent <- function (r.mat, fdr = TRUE, direct = FALSE, plot = TRUE, ...) 
{
  pcor = sm2vec(r.mat)
  indexes = sm.index(r.mat)
  colnames(indexes) = c("node1", "node2")
  w = cbind(pcor, indexes)
  if (fdr == TRUE) {
    # cat("Estimate (local) false discovery rates (partial correlations):\n")
    fdr.out = fdrtool(w[, 1], statistic = "correlation", 
                      plot = plot, ...)
    pval = fdr.out$pval
    qval = fdr.out$qval
    prob = 1 - fdr.out$lfdr
  }
  else {
    pval = rep(NA, length(w[, 1]))
    qval = pval
    prob = pval
  }
  result = cbind(w, pval, qval, prob)
  if (direct == TRUE) {
    spvar = attr(r.mat, "spv")
    if (is.null(spvar)) {
      r.mat.cor = pcor2cor(r.mat)
      spvar = 1/diag(solve(r.mat.cor))
    }
    p = length(spvar)
    r.spvar = (t(spvar %*% t(rep(1, p)))/(spvar %*% t(rep(1, 
                                                          p))))
    log.spvar = log(sm2vec(r.spvar))
    if (fdr == TRUE) {
      if (plot == TRUE) {
        dev.new()
      }
      # cat("Estimate (local) false discovery rates (log ratio of spvars):\n")
      fdr.out = fdrtool(log.spvar, statistic = "normal", 
                        plot = plot, ...)
      pval.dir = fdr.out$pval
      qval.dir = fdr.out$qval
      prob.dir = 1 - fdr.out$lfdr
    }
    else {
      pval.dir = rep(NA, length(w[, 1]))
      qval.dir = pval.dir
      prob.dir = pval.dir
    }
    result = cbind(result, log.spvar, pval.dir, qval.dir, 
                   prob.dir)
  }
  sort.idx = order(-abs(result[, 1]))
  result = as.data.frame(result[sort.idx, ])
  return(result)
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
  gn_pvalues <- network.test.edges.silent(adja_data, plot=F, verbose=F)
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

# reproduce Figure 3B (glycomics data)
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
  
  # plot
  p <- ggplot(as.data.frame(x), aes(x=cut_vec, y=-median_p))+
    geom_errorbar(aes(ymin=-upperbound, ymax=-lowerbound), colour="grey")+
    geom_line(size=1)+
    xlab("Correlation cutoff")+
    ylab("-log10(Fisher's test p-value)")+
    theme_bw() +
    ggtitle("GeneNet Cutoff optimization")
  
  # save plot to file
  pdf("Figure3B.pdf", width = 10.5, height= 10)
  print(p)
  dev.off()
  
  return(p)
  
}

# reproduce Figure 3C (glycomics data)
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

# reproduce Figure 4A (glycomics data)
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

# reproduce Figure 4B (glycomics data)
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

# reproduce Figure 4C (glycomics data)
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
  
  # plot
  p <- ggplot(as.data.frame(x), aes(x=cut_vec, y=-median_p, color=Adjacency))+
    geom_errorbar(aes(ymin=-upperbound, ymax=-lowerbound), colour="grey")+
    geom_line(size=1)+
    xlab("Correlation cutoff")+
    ylab("Fisher's test p-value (-log10)")+
    scale_colour_manual(values = c("Original"= "black","1-Sugar"="dodgerblue3", "Subclass"="seagreen3")) +
    theme_light(base_size = 15) +
    ggtitle("GeneNet Cutoff optimization")
  
  # save plot to file
  pdf("Figure4C.pdf", width = 10.5, height= 10)
  print(p)
  dev.off()
  
  return(p)
}

# function to download TCGA PANCAN12 RNA-seq data from Xena Browser
f_download_TCGA_data_from_XenaBrowser <- function(url, name_){
  # data for pancan12 mrna/clin data 
  tmp = paste0("data/pancan12_downloaded_", name_)
  # if not downloaded yet 
  if(length(intersect(list.files(), tmp))<1 ) download.file(url, paste0(getwd(),"/",tmp))
  read.table(gzfile(tmp), header = T, row.names = 1,sep = "\t")
}

# arrange TCGA PANCAN12 RNA-seq data
download_rna <- function(){
  
  # urls of data deposited in Xena Browser
  urls <- list(mrna = "https://legacy.xenahubs.net/download/TCGA.PANCAN12.sampleMap/PanCan12.3602-corrected-v3_syn1715755.gz",
               clin = "https://legacy.xenahubs.net/download/TCGA.PANCAN12.sampleMap/PANCAN12_clinicalMatrix.gz" )
  
  # download mrna data
  mrna <- t(f_download_TCGA_data_from_XenaBrowser(urls$mrna,"mrna"))
  # download clinical data 
  clin <- f_download_TCGA_data_from_XenaBrowser(urls$clin,"clin")[,c("age_at_initial_pathologic_diagnosis","gender","OS","OS.time","X_primary_disease")]
  
  # format the IDs
  rownames(mrna)=gsub("-", ".", rownames(mrna))
  rownames(clin)=gsub("-", ".", rownames(clin))
  
  # select only patients having both mrna and clin data 
  inds = sort(intersect(rownames(mrna), rownames(clin)))
  mrna = mrna[inds,]
  clin = clin[inds, ]
  
  # saniity check
  identical(rownames(mrna), rownames(clin))
  
  # create survival data
  S=Surv(time = clin$OS.time,event = clin$OS)
  
  # create final data.frame
  clinical.vars_mrna=data.frame(S=S, age=clin$age_at_initial_pathologic_diagnosis,
                                gender=clin$gender, ctype=clin$X_primary_disease)
  rownames(clinical.vars_mrna) <- rownames(clin)
  
  mrna.pan12=mrna
  
  # save data to files
  save(file = "data/pancan12.mrna",clinical.vars_mrna,mrna.pan12)
}

# knn imputation function
fimp_knn_on_subspace <- function(dta, drop.th =0.2, cr.th = 0.4, N = 10, k.min = 5, k.max = Inf, place_on = T){
  
  # missingness threshold
  yh = colSums(is.na(dta))
  D = dta
  D = as.matrix(D[, yh < nrow(D)*drop.th ])
  yh = colSums(is.na(D))
  
  imputed.values<-
    lapply( which(yh>0), function(i) {
      
      m = D[,i] # missing data
      D_ = D[,-i] # rest
      # complete candidate variables
      inds2 = colSums(is.na(D_[!is.na(m),])) == 0
      D_ = D_[, inds2]
      
      cr = as.vector(cor(m, D_,  use = "pairwise.complete.obs"))
      # constraint feature space 
      inds = abs(cr) > cr.th
      l = sum(inds)
      # cat(l)
      
      if(l < k.min){
        cr.th = abs(cr[order(abs(cr),decreasing = T)[k.min]])
        inds = abs(cr) >= cr.th
        l = sum(inds)
      }
      
      if(l > k.max){
        cr.th = abs(cr[order(abs(cr),decreasing = T)[k.max]])
        inds = abs(cr) >= cr.th
        l = sum(inds)
      }
      
      # print(l)
      stat = c(n.cand = sum(inds2), n.csig = l, cr = min(abs(cr)[inds]))
      
      D_ = D_[,inds]
      dm = pdist::pdist(D_[is.na(m), ], D_[!is.na(m), ])
      dm = t(matrix(dm@dist,dm@p,dm@n))
      
      # complete observations
      comp.obs <- which(!is.na(m))
      
      imputed<-
        apply(dm,1,function(x){
          # x= dm[1,]
          ns = order(x)[1:N]
          d = exp(-x[ns])
          ns = comp.obs[ns]
          (m[ns] %*% d)/sum(d)
          
        })
      attributes(imputed)$stats <- stat
      
      imputed
    })
  
  if(place_on){
    # missing variables
    miss.inds = which(yh>0)
    res<-
      sapply(1:length(miss.inds), function(i){
        x = D[, miss.inds[i]]
        x[is.na(x)] = imputed.values[[i]]
        x
      })
    D[,miss.inds] <- res
    return(D)
  }
  
  return(imputed.values)
  
}

# impute TCGA PANCAN12 RNA-seq data
impute_rna <- function(){
  
  # load mrna data created by download_rna()
  load("data/pancan12.mrna")
  
  cancers = levels(as.factor(clinical.vars_mrna$ctype))
  names(cancers) <- cancers
  
  # impute missing data
  mrna.pan12.imputed <- lapply(cancers, function(tag){
    print(tag)
    dta = as.matrix( mrna.pan12[clinical.vars_mrna$ctype == tag, ] )
    fimp_knn_on_subspace(dta)
  })
  # 
  # cancers = levels(clinical.vars_mrna$ctype)
  # dta = mrna.pan12[clinical.vars_mrna$ctype == cancers[1], ]
  # 
  # dta = as.matrix(dta)
  # 
  # dta_ <- fimp_knn_on_subspace(dta)
  # 
  
  # check dimensions
  t(sapply(mrna.pan12.imputed, dim))
  
  # ssanity check
  sapply(mrna.pan12.imputed, function(x) any(is.na(x)))
  
  # select clinical data
  clin.mrna <- lapply(mrna.pan12.imputed, function(x) clinical.vars_mrna[rownames(x), ])
  
  # save data to file
  save(file = "data/pancan12.mrna.imputed",  mrna.pan12.imputed, clin.mrna)
}

# correct TCGA PANCAN12 RNA-seq data for age, gender and cancer type
preprocess_rna <- function(){
  
  # load mrna data created by impute_rna()
  load("data/pancan12.mrna.imputed")
  
  # get gene names
  mrna_genes = lapply(mrna.pan12.imputed, colnames) %>% 
    unlist %>% table %>% `==`(12) %>% which %>% names
  # # string adjacency matrix
  # load("adj_string")
  # adj_string = 1*(adj_string != 0)
  # 
  # # genes in both string and tcga all pancan12  
  # common_genes = intersect(mrna_genes), colnames(adj_string))
  common_genes = mrna_genes
  
  # format clinical data
  clin.mrna = lapply(clin.mrna, function(x) {x$gender %<>% as.character; x})
  clin_mrna = do.call(rbind, clin.mrna)
  clin_mrna$gender[(clin_mrna$gender == '')] = NA
  rownames(clin_mrna) = lapply(clin.mrna, rownames) %>% unlist
  clin_mrna$gender %<>% {as.numeric(as.factor(.))}
  
  # correct data for age and gender
  mrna = 
    lapply(mrna.pan12.imputed, function(x){
      cl = clin_mrna[rownames(x), c("age","gender")] %>% as.matrix
      scale(apply(x, 2, function(x) residuals(lm(x~cl))))
    })
  mrna = mrna %>% { rn = lapply(., rownames) %>% unlist; 
  
  # select only genes common to all cancers
  x = do.call(rbind, lapply(., function(y) y[,common_genes]));
  rownames(x) = rn; x }
  rownames(mrna) = unname(rownames(mrna))
  clin = clin_mrna[rownames(mrna),]
  
  # sanity check
  identical(rownames(clin), rownames(mrna))
  
  # save data to file
  save(file = "data/pancan_mrna_imputed_corrected_combined_not_reduced_genes", mrna, clin)
}

# function to download TCGA pancan12 data from Xena Browser
f_download_STRING<- function(url, name_){
  # data for pancan12 mrna/clin data 
  tmp = paste0(name_)
  # if not downloaded yet 
  if(length(intersect(list.files(), tmp))<1 ) download.file(url, paste0(getwd(),"/",tmp))
  read.table(gzfile(tmp), header = T)
}

# download and subselect STRING PPI network
get_string_adjacency <- function(){
  
  # url for String network 
  url = "https://version-10-5.string-db.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz"
  
  # String in the form of edge list 
  ppi = f_download_STRING(url, "data/adj_stringV105")
  
  # standardize the ids
  ppi$protein1<-gsub(x=as.character(ppi$protein1),pattern = "9606.",replacement = "")
  ppi$protein2<-gsub(x=as.character(ppi$protein2),pattern = "9606.",replacement = "")
  
  # cretae adjacency matrix 
  a_ppi = get.adjacency( graph.data.frame(ppi), sparse=F)
  a_ppi[as.matrix(ppi[,1:2])]<-ppi$combined_score
  
  # adjacency matrix
  A <- as(a_ppi*0.001,"dgCMatrix")
  # save(file = "data/A", A)
  # convert ensembl protein ids to gene symbols 
  ENSPs <- unique(c(ppi$protein1,ppi$protein2))
  ENSPs <- ENSPs[!is.na(ENSPs)]
  
  set.seed(42)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mape <- getBM(filters= "ensembl_peptide_id", 
                attributes= c("ensembl_peptide_id","hgnc_symbol"),
                values= ENSPs, mart= mart)
  
  sum(mape$hgnc_symbol=="")
  
  A <- A[mape$ensembl_peptide_id, mape$ensembl_peptide_id]
  
  identical(rownames(A), mape$ensembl_peptide_id)
  identical(colnames(A), mape$ensembl_peptide_id)
  
  keep=names(which(table(mape$hgnc_symbol)==1))
  
  mape=mape[mape$hgnc_symbol%in% keep,]
  A <- A[mape$ensembl_peptide_id, mape$ensembl_peptide_id]
  colnames(A) <- rownames(A) <- mape$hgnc_symbol
  
  adj_string = A
  save(file = "data/adj_string", adj_string)
  
}

# download reactome pathways
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

# run cutoff optimization based on a list of adjacency matrices
ffcutopL <- function(adjl, data, nboots=100, cut_vec = seq(from = 0, to = 1, length = 100)){
  
  #calculate 95% confidence interval of the values (not CI of mean)
  confin_extended <- function(fisp) {
    quantile(fisp, probs = c(0.025,0.5,0.975), na.rm = T)
  }
  
  
  # cutoff optimization function
  frun_pcout <- function(data_re, cut_vec, p_ths = c(0.05,0.01), 
                         padj_methods = c("fdr","bonferroni"), return_cormat = F){
    # compute partial correlation
    cormat <- ggm.estimate.pcor(as.matrix(data_re), method = "dynamic",  verbose=FALSE)
    # compute partial correlation p-values
    gn_pvalues <- network.test.edges.silent(cormat, plot=FALSE, verbose=FALSE)
    
    # get partial correlation cutoffs for bonferroni and fdr cutoffs 
    pcor_vals<-
      padj_methods %>% lapply( function(meth) 
        p_ths %>% lapply(function(p.th)
          p.adjust(gn_pvalues$pval,meth) %>% 
            {a = gn_pvalues[. <= p.th,];a$pcor[which.max(a$pval)]} %>%
            ifelse(length(.)==0,NA,.) %>% abs
        )) %>% unlist %>% 
      { names(.) = c("fdr_0.05", "fdr_0.01", "bonf_0.05", "bonf_0.01");.}
    
    # get Fisher's p-value and statistics for each cutoff
    ps_for_cutoffs <- lapply(adjl, function(adja){
      sapply(c(pcor_vals,cut_vec), function(cutoff){
        co = as.vector(as.dist(1*(abs(cormat) >= cutoff)))
        ad = as.vector(as.dist(adja))
        fisher.test( table(factor(co,levels = 0:1),
                           factor(ad,levels = 0:1) )) %>%
          {c(p = unname(log10(.$p.value)), t = unname(.$estimate) )}
      }) %>% t
    })
    
    if(return_cormat) 
      return(list(pcor_vals = pcor_vals, pvals_for_cutoffs = ps_for_cutoffs, ppcor = cormat ))
    
    # return p-values
    list(pcor_vals = pcor_vals, pvals_for_cutoffs = ps_for_cutoffs)
  }
  
  # run calculation on original data 
  m0 = frun_pcout(data, cut_vec, return_cormat=T)
  
  # run calculation on bootstrapped data
  mboo =
    seq(nboots) %>% lapply(function(i){
      frun_pcout(data[sample(nrow(data), nrow(data), replace = TRUE),],cut_vec)
    })
  
  # to collect results, iterate over adjacency list
  rel <- rep(list(NULL), length(adjl))
  names(rel) <- names(adjl)
  for(i in seq(adjl)){
    # estimate confidence intervals
    pcutofs = sapply(mboo, function(x) x$pvals_for_cutoffs[[i]][,"p"])
    pconf_cutoffs = cbind(m0=m0$pvals_for_cutoffs[[i]][,"p"], apply( pcutofs,1, confin_extended) %>% t)
    
    stat_cutoffs = sapply(mboo, function(x) x$pcor_vals)
    # set p=1 for NAs
    # stat_cutoffs[is.na(stat_cutoffs)] = 1
    conf_stat_cutoffs = cbind(m0=m0$pcor_vals, apply( stat_cutoffs,1, confin_extended) %>% t )
    
    # name them
    colnames(pconf_cutoffs)[-1] <- colnames(conf_stat_cutoffs)[-1] <- c("L","M","U")
    rel[[i]] <-  list(ntw_cut_offs_ps = pconf_cutoffs, stat_cutoffs =conf_stat_cutoffs)
  }
  
  list(res_adjl = rel, pcormat = m0$ppcor)
}

# create adjaceny matrix from edgelist
fl2a <- function(x) ( as_adjacency_matrix(
  graph.data.frame(x,directed=FALSE),
  type="both",names=TRUE,sparse=FALSE) > 0) *1

# get adj matrijes such as s_ = S/C, or r_ = R/C
get_unique_adj <- function(adj_matrices){
  ff_ <- function(adjl, cadj){
    re = lapply(adjl, function(adj) adj %>%{.[cadj == 1] = 0;.})
    names(re) = paste0(names(adjl),"_")
    c(re,c= list(cadj))
  }
  
  ff_(adj_matrices[1:2], adj_matrices$c)
}

# get contingency table with p values, given ppcor mat, 
# reference network and, a cutoff
fctp <- function(cormat, adja, cutoff){
  co = as.vector(as.dist(1*(abs(cormat) >= cutoff)))
  ad = as.vector(as.dist(adja))
  tb = table(factor(co,levels = 0:1),
             factor(ad,levels = 0:1) )
  
  structure(tb, cx = cutoff,
            pv = log10( fisher.test(tb)$p.value ) )
}

# run cut opt procedure and gather results from all possible train, test cases
frun_and_gather_results_from_nod <- function(x, nBoo = 5, cut_vec){
  try({ 
    priors = get_unique_adj(x$A)
    # lsit of results for each network
    rxx = ffcutopL(priors , mrna[,colnames(x$A$c)], nboots = nBoo, cut_vec = cut_vec)
    
    cases = c("s_","r_","c","c")
    re<-
      cbind(tr = cases, ts = rev(cases)) %>% 
      {rownames(.) = apply(.,1, paste,collapse = ".");.} %>% apply(1, function(k){
        # priors[k[1]] #train
        # trained network 
        ntw_cut_offs_ps = rxx$res_adjl[[k[1]]]$ntw_cut_offs_ps[-(1:4),]
        tr.cutoff = cut_vec[which.min(ntw_cut_offs_ps[,"M"])]
        tr.statcutoff = rxx$res_adjl[[k[1]]]$stat_cutoffs[1,1]
        
        # same for test and train
        cormat = rxx$pcormat
        ts.adja = priors[[k[2]]] #test
        
        list( ntw = fctp(cormat, ts.adja, tr.cutoff),
              stat = fctp(cormat, ts.adja, tr.statcutoff) )
      })
    return(list(res_trts = re, res_opt = rxx))
  })
} 
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
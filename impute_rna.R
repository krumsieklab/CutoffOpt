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
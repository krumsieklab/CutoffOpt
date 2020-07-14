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
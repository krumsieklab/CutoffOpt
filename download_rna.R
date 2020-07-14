# function to download TCGA pancan12 data from Xena Browser
f_download_TCGA_data_from_XenaBrowser <- function(url, name_){
  # data for pancan12 mrna/clin data 
  tmp = paste0("data/pancan12_downloaded_", name_)
  # if not downloaded yet 
  if(length(intersect(list.files(), tmp))<1 ) download.file(url, paste0(getwd(),"/",tmp))
  read.table(gzfile(tmp), header = T, row.names = 1,sep = "\t")
}

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
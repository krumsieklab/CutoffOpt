# function to download TCGA pancan12 data from Xena Browser
f_download_STRING<- function(url, name_){
  # data for pancan12 mrna/clin data 
  tmp = paste0(name_)
  # if not downloaded yet 
  if(length(intersect(list.files(), tmp))<1 ) download.file(url, paste0(getwd(),"/",tmp))
  read.table(gzfile(tmp), header = T)
}

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
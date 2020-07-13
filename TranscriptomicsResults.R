#### Initialize ----

# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import helper functions
source("HelperFunctions.R")

# load libraries
library(magrittr)
library(parallel)
library(graphite)
library(pdist)
library(survival)
library(igraph)
library(GeneNet)
library(tidyr)
library(ggplot2)
library(scales)
library(plotly)
library(gridExtra)
library(ggrepel)
library(visNetwork)
library(dplyr)
library(ggpubr)

# check if directory "data" exists, otherwise create
if(!dir.exists("data")) {
  dir.create("data")
}

#### Load STRING adjacency ----

# download precomputed STRING adjacency (access date 06/14/2019)
load.web.file(
  url = "https://ndownloader.figshare.com/files/23821079?private_link=477d393facf01dda8355", 
  md5sum = "df8ccca5d1eecd16fe01488f490018d1", 
  outfile = "data/adj_string", 
  zipfile = F)
load("data/adj_string")

# To reproduce the precomputed STRING from scratch and load them into the script, run the lines below (running time ~1h)
# source("get_string_adjacency.R")
# library(igraph)
# library(biomaRt)
# library(Matrix)
# get_string_adjacency()
# load("data/adj_string")
# HOWEVER, please note that since packages related to gene id mapping are updated in time, depending on the access date
# the final adjacency might be slightly different than the one used in our paper. In order to reproduce our results exactly, 
# we hence recommend to load the precomputed version provided.

#### Load CORUM adjacency ----

# download precomputed CORUM adjacency (access date 06/14/2019)
load.web.file(
  url = "https://ndownloader.figshare.com/files/23821082?private_link=477d393facf01dda8355", 
  md5sum = "140a9b0f18f8417b26c5d8c61d1fe59a", 
  outfile = "data/corum", 
  zipfile = F)
load("data/corum")

# To reproduce the precomputed CORUM from scratch and load them into the script, run the lines below (running time ~1h)
# library(ganet)
# corum <- fl2a(ganet::corum)
# HOWEVER, please note that the CORUM data in the ganet package might be updated in time, and hence depending on the access date
# the adjacency might be slightly different than that used for our paper. In order to reproduce our results exactly, 
# we hence recommend to load the precomputed version provided.

#### Load TCGA PANCAN12 RNA-seq data ----

# load precomputed results?
# set to F if you want to process data and generate all results from scratch
# NOTE: the computation of all results will take several days on a regular laptop
load_precomputed <- T

# load preprocessed TCGA PANCAN12 RNA-seq data (downloaded from https://xenabrowser.net/datapages/?dataset=TCGA.PANCAN12.sampleMap%2FPanCan12.3602-corrected-v3_syn1715755&host=https%3A%2F%2Flegacy.xenahubs.net&addHub=https%3A%2F%2Flegacy.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# we additionally filtered only genes common to all 12 cancer types and and corrected for age, gender and cancer type
if(!load_precomputed){
  
  source("download_rna.R")
  source("impute_rna.R")
  source("preprocess_rna.R")
  
  # download data
  download_rna() 
  # impute data
  impute_rna()
  # preprocess data
  preprocess_rna()
  
  # load data
  load("data/pancan_mrna_imputed_corrected_combined_not_reduced_genes")
} else {
  load.web.file(
    url = "https://ndownloader.figshare.com/files/23821088?private_link=477d393facf01dda8355", 
    md5sum = "74e8e75766761a63af7165caec22152d", 
    outfile = "data/TCGA_PANCAN12_preprocessed", 
    zipfile = F)
  load("data/TCGA_PANCAN12_preprocessed")
}

#### Load Reactome pathways ----

# download precomputed Reactome pathway annotations obtained from the graphite package (access date 06/14/2019)
# the following filters were applied to subselect the pathways:
# 1. pathway contains at least 10 genes and at most 1,000 genes 
# 2. at least 50% of the genes in the pathway were found in the PANCAN12 RNA-seq data
load.web.file(
  url = "https://ndownloader.figshare.com/files/23821085?private_link=477d393facf01dda8355", 
  md5sum = "e6e631d53b24152e91c8d1c92ac1c52c", 
  outfile = "data/graphite_reactome_pathways_061419", 
  zipfile = F)
load("data/graphite_reactome_pathways_061419")

# To reproduce the precomputed reactome pathways from scratch and load them into the script, run the lines below (running time ~1h)
# source("get_reactome_pathways.R")
# get_reactome_pathways()
# load("data/graphite_reactome_pathways")
# HOWEVER, please note that these pathway annotations are updated in time, and hence depending on the graphite package updates
# they might be slightly different than those used for our paper. In order to reproduce our results exactly, 
# we hence recommend to load the precomputed version provided.

#### Set global parameters ----

# define vector of cutoffs
cut_vec = seq(from = 0, to = 1, length = 100)
# define number of bootstrapping
nboot <- 100

#### Subset adjacency matrices to only genes in data ----

# keep only part of CORUM overlapping with RNAseq data
corum = corum[intersect(colnames(corum), colnames(mrna)),
              intersect(colnames(corum), colnames(mrna))]

# keep only part of STRING overlapping with RNAseq data
adj_string = as.matrix(1*(adj_string > 0))
adj_string = adj_string[intersect(colnames(adj_string), colnames(mrna)),
                        intersect(colnames(adj_string), colnames(mrna))]

#### Create gene lists and STRING/CORUM adjacency for each REACTOME pathway ----

nods <-
  # loop over pathways
  lapply(pws, function(x){

    # extract reactome adjacency for given pathway
    radj = fl2a(x@protEdges[,c("src","dest")])
    # get overlapping genes between reactome pathway and data
    genes_r_m = intersect(colnames(mrna), colnames(radj)) %>% unique
    # restrict adjacency and data to only common genes
    radj = radj[genes_r_m,genes_r_m]
    pw_mrna = mrna[ ,genes_r_m]
    
    # subset STRING adjacency to pathway genes
    genes_r_m_s =intersect(genes_r_m, colnames(adj_string))
    hh = adj_string[genes_r_m_s, genes_r_m_s]
    sadj = radj*0
    sadj[colnames(hh),colnames(hh)] = hh
    rm(hh)
    
    # subset CORUM adjacency to pathway genes
    genes_r_m_c =intersect(genes_r_m, colnames(corum))
    hh = corum[genes_r_m_c, genes_r_m_c]
    cadj = radj*0
    cadj[colnames(hh),colnames(hh)] = hh
    rm(hh)
    
    # string - corum, reactome - corum
    s_c_adj = sadj %>%{.[cadj == 1] = 0;.}
    r_c_adj = radj %>%{.[cadj == 1] = 0;.}
    
    # function to calculate mean density
    fpwden <- function(A) c(dn =sum(A)/((nrow(A)^2)-nrow(A)))
    # function to calculate mean degree
    fpwdeg <- function(A) c(dg = (sum(A)/2)/(nrow(A)))
    
    # densities
    adj_list = list(s = sadj, r = radj, c = cadj, s_ =s_c_adj, r_ = r_c_adj)
    list( A =  list(s = sadj, r = radj, c = cadj), #list(s_ = s_c_adj, r_ = r_c_adj, c = cadj),
          dn = c( adj_list %>% sapply(fpwden),
                  n = nrow(cadj),
                  adj_list %>% sapply(fpwdeg))
    )
  })

#### Run simulation / Load precomputed data ----

if(load_precomputed) {
  
  # load precalculated data (nboot=100)
  load.web.file(
    url = "https://ndownloader.figshare.com/files/23821091?private_link=477d393facf01dda8355", 
    md5sum = "3e494b6d0a5bdc777d955882bb72ad8a", 
    outfile = "data/TranscriptomicsCutoptResults", 
    zipfile = F)
  load("data/TranscriptomicsCutoptResults")
  ress <- ress_100
  rm(ress_100)
  
} else {
  
  tic()
  # loop over pathways
  ress <- lapply(nods, function(x) {
    # run cutoff optimization on STRING and compute overlap on CORUM
    frun_and_gather_results_from_nod(x=x,nBoo = nboot, cut_vec = cut_vec)
  })
  toc()
  
}

#### Figure 6 ----

# only consider pathways < 1000 genes
q <- lapply(nods, function(x) x$dn[which(names(x$dn)=="n")]) %>% unlist 
del.pw <- which(q > 1000)

ress[del.pw] <- NULL
nods[del.pw] <- NULL

# create data for plot
overview_plot <- list()
overview_plot$name <- names(ress)
overview_plot$dim <- lapply(nods, function(x) x$dn[which(names(x$dn)=="n")]) %>% unlist

# extract Fisher's test p-value at optimal cutoff using STRING
overview_plot$p_s_opt <- lapply(ress, function(x) {min(x$res_opt$res_adjl$s_$ntw_cut_offs_ps[-c(1:4),"m0"][which.min(x$res_opt$res_adjl$s_$ntw_cut_offs_ps[-c(1:4),"M"])])}) %>% unlist
# extract Fisher's test p-value at FDR cutoff
overview_plot$p_s_fdr <- lapply(ress, function(x) x$res_opt$res_adjl$s_$ntw_cut_offs_ps[1,1]) %>% unlist
# extract Fisher's test p-value using optimal cutoff tested on CORUM
overview_plot$p_s_c_opt <- lapply(ress, function(x) sapply(x$res_trts$s_.c, attr, "pv")[1]) %>% unlist
# extract Fisher's test p-value using FDR cutoff tested on CORUM
overview_plot$p_s_c_fdr <- lapply(ress, function(x) sapply(x$res_trts$s_.c, attr, "pv")[2]) %>% unlist

overview_plot$or_s_c_opt <- lapply(ress, function(x) fisher.test(x$res_trts$s_.c$ntw)$estimate) %>% unlist
overview_plot$or_s_c_fdr <- lapply(ress, function(x) fisher.test(x$res_trts$s_.c$stat)$estimate) %>% unlist

overview_plot <- as.data.frame(overview_plot)

# create summary plot of optimization results on STRING
inds1 <- overview_plot$p_s_opt < log10(0.01/dim(overview_plot)[1])
hh1 <- ggplot(overview_plot[!inds1,], aes(x=-p_s_opt, y=-p_s_fdr, size=dim, label=name)) + 
  geom_point(fill = "gray", alpha=0.3, pch=21, color="gray") +
  geom_point(data=overview_plot[inds1,], fill="lightsteelblue2",
             aes(size=dim), alpha=0.7, pch=21, color= "black") +
  scale_x_continuous(trans="sqrt", breaks=pretty_breaks(10), expand=c(0,0), limits=c(0,10+max(c(-overview_plot$p_s_opt,-overview_plot$p_s_fdr)))) +
  scale_y_continuous(trans="sqrt", breaks=pretty_breaks(10), expand=c(0,0), limits=c(0,10+max(c(-overview_plot$p_s_opt,-overview_plot$p_s_fdr)))) +
  geom_text_repel(data=overview_plot[overview_plot$name=="Axon guidance",],
                  size=5) +
  geom_hline(data = NULL, yintercept=-log10(0.01/dim(overview_plot)[1]), colour='grey50', linetype=2,
             na.rm = FALSE, show.legend = NA)+
  geom_vline(data = NULL, xintercept=-log10(0.01/dim(overview_plot)[1]), colour='grey50', linetype=2,
             na.rm = FALSE, show.legend = NA)+
  geom_abline(data = NULL, slope=1, intercept=0, colour='red', linetype=2,
              na.rm = FALSE, show.legend = NA)+
  theme_minimal() +
  labs(size="Pathway size", fill="Pathway size") +
  ggtitle("Overlap with STRING reference")+ 
  xlab("-log10(Fisher's p-value Optimal Cutoff based on STRING)")+
  ylab("-log10(Fisher's p-value FDR 0.05 Cutoff)")

# create summary plot of validation results on CORUM
inds2 <- overview_plot$p_s_c_opt < log10(0.01/dim(overview_plot)[1])
hh2 <- ggplot(overview_plot[inds1 & !inds2,], aes(x=-p_s_c_opt, y=-p_s_c_fdr, size=dim, label=name)) + 
  geom_point(fill = "gray", alpha=0.3, pch=21, color="gray") +
  geom_point(data=overview_plot[inds1 & inds2,], fill="lightsteelblue2",
             aes(size=dim), alpha=0.7, pch=21, color= "black") +
  scale_x_continuous(trans="sqrt", breaks=pretty_breaks(10), expand=c(0,0), limits=c(0,10+max(c(-overview_plot$p_s_c_opt,-overview_plot$p_s_c_fdr)))) +
  scale_y_continuous(trans="sqrt", breaks=pretty_breaks(10), expand=c(0,0), limits=c(0,10+max(c(-overview_plot$p_s_c_opt,-overview_plot$p_s_c_fdr)))) +
  geom_text_repel(data=overview_plot[overview_plot$name=="Axon guidance",],
                  size=5) +
  geom_hline(data = NULL, yintercept=-log10(0.01/dim(overview_plot)[1]), colour='grey50', linetype=2,
             na.rm = FALSE, show.legend = NA)+
  geom_vline(data = NULL, xintercept=-log10(0.01/dim(overview_plot)[1]), colour='grey50', linetype=2,
             na.rm = FALSE, show.legend = NA)+
  geom_abline(data = NULL, slope=1, intercept=0, colour='red', linetype=2,
              na.rm = FALSE, show.legend = NA)+
  theme_minimal() +
  labs(size="Pathway size", colour="Optimization P-value") +
  ggtitle("Overlap with CORUM reference")+ 
  xlab("-log10(Fisher's p-value Optimal Cutoff based on STRING)")+
  ylab("-log10(Fisher's p-value FDR 0.05 Cutoff)")

# save plot to file
op <- ggarrange(hh1, hh2, ncol=2, common.legend = T, legend="right", labels = "AUTO")
ggsave("Figure6.pdf", op, width = 10, height = 5)

# save interactive html version
op_int <- subplot(ggplotly(hh1, tooltip = c("label")),
                  ggplotly(hh2, tooltip = c("label")))
op_int <- op_int %>% layout(annotations = list(
  list(x = 0.03, y = 1.03, text = "Overview optimization on String", showarrow = F, xref='paper', yref='paper'),
  list(x = 0.8, y = 1.03, text = "Overview validation on CORUM (opt on String)", showarrow = F, xref='paper', yref='paper')),
  title="")
htmlwidgets::saveWidget(as_widget(op_int), "SupplementaryData1.html")

#### Supplementary Figure S8 ----

# collect all optimization results for STRING
cv <- rep(cut_vec, times = length(ress))

# get pathway size
pw.dim <- lapply(nods, function(x) x$dn[which(names(x$dn)=="n")])

# create label for each pathway
nn <- lapply(names(ress), function(x) paste(x," (",pw.dim[[x]],")", sep = "")) %>% unlist
nn <- rep(nn, each = length(cut_vec))

# get mean Fisher's test p-value for the statistical cutoffs across the bootstrapping
pval.opt <- lapply(ress, function(x) min(x$res_opt$res_adjl$s_$ntw_cut_offs_ps[,"M"])) %>% unlist
# assess if the p-value is significant
is.sign <- pval.opt < log10(0.01/length(ress))
is.sign <- rep(is.sign, each = length(cut_vec))

oo <- lapply(ress, function(x) x$res_opt$res_adjl$s$ntw_cut_offs_ps[-c(1:4,105),]) 
oo <- do.call("rbind", oo)

data_plot <- cbind.data.frame(name=nn, cut_vec=cv, sign=is.sign, oo)

# collect all statistical cutoffs
cut_stat <- lapply(ress, function(x) x$res_opt$res_adjl$s_$stat_cutoffs[,1])
cut_stat <- as.data.frame(do.call("rbind", cut_stat))
cut_stat <- stack(cut_stat)[1]

cut_stat$test <- rep(c("fdr_0.05","fdr_0.01","bonf_0.05","bonf_0.01"), each = length(ress))

cut_stat$name <- rep(lapply(names(ress), function(x) paste(x," (",pw.dim[[x]],")", sep = "")) %>% unlist, times=4)

g <- ggplot(data_plot, aes(x=cut_vec, y=-M, color=sign))+
  geom_errorbar(aes(ymin=-U, ymax=-L,color=sign))+
  scale_colour_manual(values = c("grey", "seagreen3"))+
  geom_line(size=1) +
  xlab("Partial correlation cutoff")+
  ylab("-log10(Fisher's p-value)")+
  theme_light(base_size=10) +
  theme(strip.text = element_text(size=6))+
  theme_minimal()+
  geom_vline(data=cut_stat,aes(xintercept=values, linetype=test))+
  scale_linetype_manual(values=c(4,3,2,1),
                        name="Statistical cutoffs",
                        breaks=c("fdr_0.05","fdr_0.01","bonf_0.05","bonf_0.01"),
                        labels=c("FDR 0.05", "FDR 0.01", "Bonf 0.05", "Bonf 0.01")) +
  geom_hline(yintercept = -log10(0.01/length(ress)), color ="red", linetype="dashed",size=0.25)+
  facet_wrap(name~., ncol = 8, scales = "free_y")

ggsave(filename='SupplementaryFigure8.pdf', plot=g, width=30, height=40, units='in', limitsize=FALSE)

#### Figure 7A ----

# create optimization plots for each significant pathway separately
g_single <- lapply(names(ress)[pval.opt < log10(0.01/length(ress))] %>% {names(.)=.;.}, function(x) {
  
  # collect optimization results
  oo <- ress[[x]]$res_opt$res_adjl$s$ntw_cut_offs_ps[-c(1:4),]
  data_plot <- as.data.frame(cbind(cut_vec=cut_vec, oo))
  # collect all statistical cutoffs
  cut_stat <- as.data.frame(ress[[x]]$res_opt$res_adjl$s_$stat_cutoffs)
  colnames(cut_stat) <- "values"
  
  g <- ggplot(data_plot, aes(x=cut_vec, y=-M))+
    geom_errorbar(aes(ymin=-U, ymax=-L),color="grey50")+
    geom_line(size=1) +
    xlab("Partial correlation cutoff")+
    ylab("-log10(Fisher's p-value)")+
    theme_light(base_size=10) +
    theme(strip.text = element_text(size=6))+
    geom_vline(aes(xintercept=m0, linetype=rownames(ress[[x]]$res_opt$res_adjl$c$stat_cutoffs)), 
               data=as.data.frame(ress[[x]]$res_opt$res_adjl$c$stat_cutoffs))+
    scale_linetype_manual(values=c(3,4,2,1),
                          name="Statistical cutoffs",
                          breaks=c("fdr_0.05","fdr_0.01","bonf_0.05","bonf_0.01"),
                          labels=c("FDR 0.05", "FDR 0.01", "Bonf 0.05", "Bonf 0.01")) +
    geom_hline(yintercept = -log10(0.01/length(ress)), color ="red", linetype="dashed",size=0.25)+
    ggtitle(x)
  
  list(opt_plot=g)
  
})
ggsave(filename="Figure7A.pdf", plot=g_single$`Axon guidance`$opt_plot, width=10, height=7, units='in')


#### Generate networks on STRING optimization ----

# save results for interactive .Rmd (TranscriptomicsNetworks.Rmd)
data <- lapply(ress, function(x) {
  x$res_trts$r_.c <- NULL
  x$res_trts$c.r_ <- NULL
  x$res_trts$c.s_ <- NULL
  
  x$res_opt$res_adjl$r_ <- NULL
  x$res_opt$res_adjl$c <- NULL
  x
})

# select significant pathways
pws <- ress[names(ress)[pval.opt < log10(0.01/length(ress))]]

net_data <- 
  # loop over the significant pathways
  lapply(names(pws), function(x) {
    
    # adjacency matrices
    priors <- get_unique_adj(nods[[x]]$A)
    
    # optimal network (STRING optimization)
    net_opt <- (abs(ress[[x]]$res_opt$pcormat) > cut_vec[which.min(ress[[x]]$res_opt$res_adjl$s$ntw_cut_offs_ps[-c(1:4),"M"])])*1
    diag(net_opt) <- 0
    # color edges according to adjacency matrices
    net_opt_edges_col <- net_opt*0
    net_opt_edges_col[net_opt == 1 & priors$s_==1] <- 1
    net_opt_edges_col[net_opt == 1 & priors$s_==0 & priors$c==1] <- 2
    net_opt_edges_col[net_opt == 1 & priors$s_==0 & priors$c==0] <- 3
    
    # FDR 0.05 network
    net_fdr <- (abs(ress[[x]]$res_opt$pcormat) > ress[[x]]$res_opt$res_adjl$s$stat_cutoffs[1,1])*1
    diag(net_fdr) <- 0
    # color edges according to adjacency matrices
    net_fdr_edges_col <- net_fdr*0
    net_fdr_edges_col[net_fdr == 1 & priors$s_==1] <- 1
    net_fdr_edges_col[net_fdr == 1 & priors$s_==0 & priors$c==1] <- 2
    net_fdr_edges_col[net_fdr == 1 & priors$s_==0 & priors$c==0] <- 3
    
    # STRING adjacency network
    g_s <- graph.adjacency(priors$s_, mode = "undirected", weighted = NULL,
                           diag = FALSE)
    ee_s <- get.data.frame(graph.adjacency(priors$s_,weighted=TRUE, mode = "undirected"))
    
    ee_s$width <- 3
    ee_s$color <- "black"
    nn_s <- as.data.frame(cbind(id=colnames(priors$s), label = colnames(priors$s), title = colnames(priors$s)))
    
    # CORUM adjacency network
    g_c <- graph.adjacency(priors$s_, mode = "undirected", weighted = NULL,
                           diag = FALSE)
    ee_c <- get.data.frame(graph.adjacency(priors$c,weighted=TRUE, mode = "undirected"))
    
    ee_c$width <- 3
    ee_c$color <- "green"
    nn_c <- as.data.frame(cbind(id=colnames(priors$c), label = colnames(priors$c), title = colnames(priors$c)))
    
    # opacity parameter
    a <- 0.5
    
    g_opt <- graph.adjacency(net_opt_edges_col, mode = "undirected", weighted = TRUE,
                             diag = FALSE)
    ee_opt <- get.data.frame(graph.adjacency(net_opt_edges_col,weighted=TRUE, mode = "undirected"))
    ee_opt$width <- 3
    ee_opt$color.color[ee_opt$weight==3] <- "#C0C0C0" # the opacity parameter requires html colors
    ee_opt$color.opacity <- 1
    ee_opt$color.opacity[ee_opt$weight==3] <- a
    ee_opt$color[ee_opt$weight==2] <- "green"
    ee_opt$color[ee_opt$weight==1] <- "black"
    ee_opt$value <- 3 # E(g_fdr)$weight
    nn_opt <- as.data.frame(cbind(id=colnames(net_opt_edges_col), label = colnames(net_opt_edges_col), title = colnames(net_opt_edges_col)))
    
    
    g_fdr <- graph.adjacency(net_fdr_edges_col, mode = "undirected", weighted = TRUE,
                             diag = FALSE)
    ee_fdr <- get.data.frame(graph.adjacency(net_fdr_edges_col,weighted=TRUE, mode = "undirected"))
    ee_fdr$width <- 3
    ee_fdr$color.color[ee_fdr$weight==3] <- "#C0C0C0" # the opacity parameter requires html colors
    ee_fdr$color.opacity <- 1
    ee_fdr$color.opacity[ee_fdr$weight==3] <- a
    ee_fdr$color[ee_fdr$weight==2] <- "green"
    ee_fdr$color[ee_fdr$weight==1] <- "black"
    ee_fdr$value <- 3 # E(g_fdr)$weight
    nn_fdr <- as.data.frame(cbind(id=colnames(net_fdr_edges_col), label = colnames(net_fdr_edges_col), title = colnames(net_fdr_edges_col)))
    
    
    #Get the coordinates of the Nodes from the sparsest network
    Coords <- layout_with_fr(g_fdr) %>% 
      as_tibble %>%
      bind_cols(data_frame(names = names(V(g_fdr))))
    
    #Set the coordinates of other networks
    NetCoords_opt <- data_frame(names = names(V(g_opt))) %>%
      left_join(Coords, by= "names")
    
    NetCoords_s <- data_frame(names = names(V(g_s))) %>%
      left_join(Coords, by= "names")
    
    NetCoords_c <- data_frame(names = names(V(g_c))) %>%
      left_join(Coords, by= "names")
    
    list(g_fdr=g_fdr, ee_fdr=ee_fdr, nn_fdr=nn_fdr, Coords=Coords, 
         g_opt=g_opt, ee_opt=ee_opt, nn_opt=nn_opt, NetCoords_opt=NetCoords_opt, 
         g_s=g_s, ee_s=ee_s, nn_s=nn_s, NetCoords_s=NetCoords_s, 
         g_c=g_c, ee_c=ee_c, nn_c=nn_c, NetCoords_c=NetCoords_c,
         cut_opt = cut_vec[which.min(ress[[x]]$res_opt$res_adjl$s_$ntw_cut_offs_ps[-c(1:4),"m0"])],
         pval_opt_train = min(ress[[x]]$res_opt$res_adjl$s_$ntw_cut_offs_ps[-c(1:4),"m0"]),
         pval_opt_val = sapply(ress[[x]]$res_trts$s_.c, attr, "pv")[1],
         cut_fdr05 = ress[[x]]$res_opt$res_adjl$s_$stat_cutoffs[1,1],
         pval_fdr05_train = ress[[x]]$res_opt$res_adjl$s_$ntw_cut_offs_ps[1,1],
         pval_fdr05_val = sapply(ress[[x]]$res_trts$s_.c, attr, "pv")[2]
    )
  })
names(net_data) <- names(pws)

save(data, net_data, file = "data/TranscriptomicsNetworkData.Rds")

#### Figure 7B-E ----

pws <- "Axon guidance"

# save optimal network
visNetwork(nodes=net_data[[pws]]$nn_opt, edges = net_data[[pws]]$ee_opt) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = as.matrix(net_data[[pws]]$NetCoords_opt[,2:3])) %>%
  visSave(file="Figure7B.html", selfcontained = TRUE, background = "white")

# save FDR network
visNetwork(nodes=net_data[[pws]]$nn_fdr, edges = net_data[[pws]]$ee_fdr) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = as.matrix(net_data[[pws]]$Coords[,1:2])) %>%
  visSave(file="Figure7C.html", selfcontained = TRUE, background = "white")

# save STRING network
visNetwork(nodes=net_data[[pws]]$nn_s, edges = net_data[[pws]]$ee_s) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = as.matrix(net_data[[pws]]$NetCoords_s[,2:3])) %>%
  visSave(file="Figure7D.html", selfcontained = TRUE, background = "white")

# save CORUM network
visNetwork(nodes=net_data[[pws]]$nn_c, edges = net_data[[pws]]$ee_c) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = as.matrix(net_data[[pws]]$NetCoords_c[,2:3])) %>%
  visSave(file="Figure7E.html", selfcontained = TRUE, background = "white") 
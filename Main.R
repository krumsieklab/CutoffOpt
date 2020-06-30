#### Initialize ----

# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import functions
source("HelperFunctions.R")
source("Figure3B.R")
source("Figure3C.R")
source("Figure4A.R")
source("Figure4B.R")
source("Figure4C.R")

# load libraries
library(ggplot2)
library(GeneNet)
library(gplots)
library(colorRamps)
library(tictoc)

#### Load Data ----

# download preprocessed glycomics data
file <- "SupplementalMaterial_DatasetS1_PreprocessedData.xls"
load.web.file(
  url="https://ndownloader.figshare.com/files/9160348",
  md5sum = '32840138512bc5a257212677ec3dc3c8',
  outfile = file
)

# load preprocessed Korčula2013 data, corrected for age and gender
df <- read_excel(file, sheet = "Korčula2013_residuals", col_names = T) %>% as.data.frame()
data <- df[,2:dim(df)[2]]
rownames(data) <- df[,1]

# load prior knowledge adjacency matrix
adja <- read.csv("data/InferredPathway.csv", header = TRUE, sep = ";", dec = ",", row.names = 1)
adja[is.na(adja)] <- 0

#### Set global parameters ----

# cutoff vector
cut_vec <- seq(from = 0, to = 1, length = 100)
# number of bootstrapping
nboot <-10

#### Figure 3B: optimization curve ----

tic()
# a file called "Figure3B.pdf" will be created in the wd
Figure3B(cut_vec=cut_vec, data=data, adja=adja, nboot=nboot)
toc()

#### Figure 3C: cutoff vs. sample size heatmap ----

# create vector of sample sizes
size_step=10
data_sizes <- seq(from = 10, to = nrow(data), by = size_step)
if(nrow(data)%%size_step != 0){
  datasizes <- c(data_sizes, nrow(data))
}

tic()
# a file called "Figure3C.pdf" will be created in the wd
Figure3C(cut_vec=cut_vec, data=data, adja=adja, nboot=nboot, data_sizes=data_sizes)
toc()

#### Figure 4A: simulated partial prior knowledge ----

# create vector of percentages
percentages <- seq(from = 0, to = 0.9, length = 10)

tic()
# a file called "Figure4A.pdf" will be created in the wd
Figure4A(cut_vec=cut_vec, data=data, adja=adja, nboot=nboot, percentages=percentages)
toc()

#### Figure 4B: simulated incorrect prior knowledge ----

# create a vector of number of edge swaps 
nswap <- c(0,1:10,seq(15,50,5))

tic()
# a file called "Figure4B.pdf" will be created in the wd
Figure4B(cut_vec, data, adja, nboot, nswap)
toc()

#### Figure 4C: coarse prior knowledge ----

# create block adjacency
adja_block <- matrix(0L, nrow = dim(adja)[1], ncol = dim(adja)[2]) 
adja_block[1:20, 1:20] <- 1
adja_block[21:40, 21:40] <- 1
adja_block[41:50, 41:50] <- 1

# load 1-sugar-addition adjacency
adja_1s <- read.csv("data/adja_1sugar.csv",header = TRUE,sep = ";",dec = ",",row.names = 1)
adja_1s[is.na(adja_1s)] <- 0

tic()
# a file called "Figure4C.pdf" will be created in the wd
Figure4C(cut_vec, data, adja, adja_block, adja_1s, nboot)
toc()
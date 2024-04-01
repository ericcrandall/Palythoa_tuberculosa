# usage: Rscript run_ConStruct_K.R K T n.chains n.iter
# where K is the number of layers (i.e. Structure clusters) to be included in the analysis
# and T (or F) indicates whether to run the spatial model
# n.chains is the number of MCMC chains

library(tidyverse)
library(conStruct)

args = commandArgs(trailingOnly=TRUE)


if (length(args) < 4) {
  stop("Must supply K, spatial(T/F), n.chains, n.iter in this order", call.=FALSE)
}
#putting a . in front of an object exempts it from the rm(ls()) statements
K <- args[1] %>% as.numeric()
spatial <- args[2] %>% as.logical()
n.chains <- args[3] %>% as.numeric()
n.iter <- args[4] %>% as.numeric()

if(spatial){
  prefix <- "sp"
}else{
  prefix <- "nonsp"
}

## Read and Convert Data
print("Reading In The Data and Converting it for use by ConStruct")
load("Ptub_ConStruct.Rdata") 
popmap <- read_tsv("popmap.tsv")

## re-arrange the names in the infile to match the order of the names in popmap
freqs <- freqs[match(popmap$IndID,str_remove(row.names(freqs),"_1")),]

localities <- read_csv("hawaii_vertices.csv") %>% 
  filter(sampled == "black") %>% map_df(rev)
localities$island[3] <- "Kamo"
localities$island[6] <- "O'ahu"

popmap <- popmap %>% left_join(localities, by = join_by(Name == island ))

latlongs <- matrix(c(popmap$long, popmap$lat), ncol = 2)

#load pairwise genetic distances between individuals
ind_geog_dists <- read.csv(file = "ptub_ind_geog_matrix.csv",header =T) %>% as.matrix()



## Run a spatial model
print(c("Running a conStruct model for K = ",K, "spatial = ", spatial))
spatial9 <- conStruct(spatial = T, K = K,
                      freqs = freqs,
                      geoDist = ind_geog_dists,
                      coords = latlongs,
                      prefix = paste(prefix,"K",K, sep = "_"),
                      n.iter = n.iter,
                      n.chains = n.chains)


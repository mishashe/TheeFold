#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
sp3 <- args[3]

# sp1="Escherichia_coli"
# sp2="Klebsiella_pneumoniae"
# sp3="Salmonella_enterica"

library(stringr)
library(foreach)
library(iterators)
library(doParallel)
require(intervals)
library(sets)
require(data.table)
#library(intervalaverage)
library(dplyr)
library(plotrix)
library(iterators)
library(rhdf5)
library(seqinr);
library(questionr)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp(paste0("/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/ThreeFold/src/intersect_weighted.cpp"))

registerDoParallel(128)
combine_tables <- function(table1,table2)
{
  names12 <- unique(c(names(table1),names(table2)))
  table1 <- table1[names12]
  setdiff21 <- setdiff(names(table2),names(table1))
  if (length(setdiff21)>0)
  {
    names(table1)[is.na(names(table1))] <- setdiff21
    table1[setdiff21] <- 0
  }
  table2 <- table2[names12]
  setdiff12 <- setdiff(names(table1),names(table2))
  if (length(setdiff12)>0)
  {
    names(table2)[is.na(names(table2))] <- setdiff12
    table2[setdiff12] <- 0
  }
  return(table1[names12] + table2[names12])
}

# dir <- "/home/misha/Documents/Development/ThreeFold/" # Bechet
dir <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/"

system(paste0("mkdir -p ",dir,"plots/"))

# sp1 <- "Escherichia_coli"
# sp2 <- "Klebsiella_pneumoniae"
# sp3 <- "Serratia"

system(paste0("mkdir -p ", dir,"data/processed/",sp3,"_",sp2,"_",sp1,"/"))
group1 <- h5ls(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"))$group[2]
group2 <- h5ls(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"))$group[2]

headers1 <- h5ls(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"))$name; headers <- unique(headers1)
headers2 <- h5ls(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers2))
headers <- headers[!headers %in% c("comp1","comp2","comp")]
rTable <- foreach(header=headers, .inorder=FALSE, .combine=combine_tables) %dopar%
{
  intervals1 <- h5read(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"), paste0(group1,"/",header))
  intervals2 <- h5read(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"), paste0(group2,"/",header))
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  wtd.table(as.integer(intersects$V2-intersects$V1+1),w=as.numeric(intersects$w))
}
write.table(data.frame(r=as.integer(names(rTable)),m=as.numeric(rTable)),paste0(dir,"data/processed/",sp3,"_",sp2,"_",sp1,"/",sp3,"_",sp2,"_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)






#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
sp3 <- args[3]
sp4 <- args[4]
sp5 <- args[5]
sp6 <- args[6]
sp7 <- args[7]
sp8 <- args[8]
sp9 <- args[9]
#species=c("Cronobacter", "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_hormaechei", "Serratia", "Salmonella_enterica", "Citrobacter", "Raoultella", "Vibrio", "Lelliottia")
#sp1 <- species[1]; sp2 <- species[2]; sp3 <- species[3]; sp4 <- species[4]; sp5 <- species[5]; sp6 <- species[6]; sp7 <- species[7]; sp8 <- species[8]; sp9 <- species[9]; 
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
Rcpp::sourceCpp(paste0("~/HGTnew/multi_comparisons/ThreeFold/src/intersect_weighted.cpp"))
dir = "~/HGTnew/multi_comparisons/"



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

system(paste0("mkdir -p ",dir,"plots/"))


system(paste0("mkdir -p ",dir,"data/processed/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"/"))
group1 <- h5ls(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"))$group[2]
group2 <- h5ls(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"))$group[2]
group3 <- h5ls(paste0(dir,"data/processed/",sp4,"_",sp1,"/",sp4,"_",sp1,".h5"))$group[2]
group4 <- h5ls(paste0(dir,"data/processed/",sp5,"_",sp1,"/",sp5,"_",sp1,".h5"))$group[2]
group5 <- h5ls(paste0(dir,"data/processed/",sp6,"_",sp1,"/",sp6,"_",sp1,".h5"))$group[2]
group6 <- h5ls(paste0(dir,"data/processed/",sp7,"_",sp1,"/",sp7,"_",sp1,".h5"))$group[2]
group7 <- h5ls(paste0(dir,"data/processed/",sp8,"_",sp1,"/",sp8,"_",sp1,".h5"))$group[2]
group8 <- h5ls(paste0(dir,"data/processed/",sp9,"_",sp1,"/",sp9,"_",sp1,".h5"))$group[2]



headers1 <- h5ls(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"))$name; headers <- unique(headers1)
headers2 <- h5ls(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers2))
headers3 <- h5ls(paste0(dir,"data/processed/",sp4,"_",sp1,"/",sp4,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers3))
headers4 <- h5ls(paste0(dir,"data/processed/",sp5,"_",sp1,"/",sp5,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers4))
headers5 <- h5ls(paste0(dir,"data/processed/",sp6,"_",sp1,"/",sp6,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers5))
headers6 <- h5ls(paste0(dir,"data/processed/",sp7,"_",sp1,"/",sp7,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers6))
headers7 <- h5ls(paste0(dir,"data/processed/",sp8,"_",sp1,"/",sp8,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers7))
headers8 <- h5ls(paste0(dir,"data/processed/",sp9,"_",sp1,"/",sp9,"_",sp1,".h5"))$name; headers <- unique(intersect(headers,headers8))

headers <- unique(headers[!headers %in% c("comp1","comp2","comp")])
rTable <- foreach(header=headers, .inorder=FALSE, .combine=combine_tables) %dopar%
{
  intervals1 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"), paste0(group1,"/",header))))
  intervals2 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"), paste0(group2,"/",header))))
  intervals3 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp4,"_",sp1,"/",sp4,"_",sp1,".h5"), paste0(group3,"/",header))))
  intervals4 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp5,"_",sp1,"/",sp5,"_",sp1,".h5"), paste0(group4,"/",header))))
  intervals5 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp6,"_",sp1,"/",sp6,"_",sp1,".h5"), paste0(group5,"/",header))))
  intervals6 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp7,"_",sp1,"/",sp7,"_",sp1,".h5"), paste0(group6,"/",header))))
  intervals7 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp8,"_",sp1,"/",sp8,"_",sp1,".h5"), paste0(group7,"/",header))))
  intervals8 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp9,"_",sp1,"/",sp9,"_",sp1,".h5"), paste0(group8,"/",header))))
  
  
  intervals1 <- (intervals1 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals2 <- (intervals2 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals3 <- (intervals3 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals4 <- (intervals4 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals5 <- (intervals5 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals6 <- (intervals6 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals7 <- (intervals7 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals8 <- (intervals8 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  
  
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals3)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals4)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals5)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals6)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals7)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals8)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  wtd.table(as.integer(intersects$V2-intersects$V1+1),w=as.numeric(intersects$w))
}
write.table(data.frame(r=as.integer(names(rTable)),m=as.numeric(rTable)),paste0(dir,"data/processed/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)


# ONLY TO SAVE INTERSECTS
registerDoParallel(20)
intersects <- foreach(header=headers, .inorder=FALSE, .combine=rbind) %dopar%
{
  intervals1 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"), paste0(group1,"/",header))))
  intervals2 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp3,"_",sp1,"/",sp3,"_",sp1,".h5"), paste0(group2,"/",header))))
  intervals3 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp4,"_",sp1,"/",sp4,"_",sp1,".h5"), paste0(group3,"/",header))))
  intervals4 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp5,"_",sp1,"/",sp5,"_",sp1,".h5"), paste0(group4,"/",header))))
  intervals5 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp6,"_",sp1,"/",sp6,"_",sp1,".h5"), paste0(group5,"/",header))))
  intervals6 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp7,"_",sp1,"/",sp7,"_",sp1,".h5"), paste0(group6,"/",header))))
  intervals7 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp8,"_",sp1,"/",sp8,"_",sp1,".h5"), paste0(group7,"/",header))))
  intervals8 <- as.data.table(as.matrix(h5read(paste0(dir,"data/processed/",sp9,"_",sp1,"/",sp9,"_",sp1,".h5"), paste0(group8,"/",header))))
  
  
  intervals1 <- (intervals1 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals2 <- (intervals2 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals3 <- (intervals3 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals4 <- (intervals4 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals5 <- (intervals5 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals6 <- (intervals6 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals7 <- (intervals7 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  intervals8 <- (intervals8 %>% group_by_all() %>% summarise(COUNT = n(), .groups="keep"))
  
  
  intersects <- as.data.table(intersect_weighted(as.matrix(intervals1), as.matrix(intervals2)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals3)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals4)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals5)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals6)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals7)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects <- as.data.table(intersect_weighted(as.matrix(intersects), as.matrix(intervals8)))
  intersects <- (intersects %>% group_by_at(vars(V1,V2)) %>% summarise(w = sum(V3), .groups="keep"))
  intersects$header <- header
  intersects$length <- as.integer(intersects$V2-intersects$V1+1)
  intersects
}
intersects <- data.frame(header=intersects$header,start=intersects$V1,stop=intersects$V2,w=intersects$w,length=intersects$length)
intersects <- intersects[intersects$length>=1000,]
intersects <- intersects[with(intersects, order(length, -w)),]
write.table(intersects,paste0(dir,"data/processed/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"_intersects.bed"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

# EXTRACT SEQUENCES
library(seqinr)
seqs <- read.fasta("/scratch/ws1/msheinman-msheinman/ThreeFold/data/external/Cronobacter.fna")
locs <- read.delim(paste0(dir,"data/processed/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"_intersects.bed"), header=F, sep="\t") # Tabs? Put sep="\t"
outf <- file(paste0(dir,"data/processed/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"/",sp9,"_",sp8,"_",sp7,"_",sp6,"_",sp5,"_",sp4,"_",sp3,"_",sp2,"_",sp1,"_intersects.fa"), 'w')

doitall <- function(x) {
  seq_id <- x[[1]]; start <- x[[2]]; stop <- x[[3]]; seq <- seqs[[seq_id]]
  seqv   <- seq[start:stop]
  header <- paste(sep="", ">", attr(seq, "name"), "_", start, "_", stop, "\n")
  cat(file=outf, header)
  cat(file=outf, toupper(paste(sep="", collapse="", seqv)), "\n")
}
apply(locs, 1, doitall)
close(outf)




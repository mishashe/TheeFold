#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]
#rTh <- args[3]  commented by FM
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
library(tidyr)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp(paste0("~/HGTnew/multi_comparisons/ThreeFold/src/intersect_weighted.cpp"))



# dir <- "/home/misha/Documents/Development/ThreeFold/"
#dir <- "/scratch/ws1/msheinman-msheinman/ThreeFold/"
dir <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/"
# sp1 <- "Escherichia_coli"
# sp2 <- "Klebsiella_pneumoniae"
# sp2 <- "Salmonella_enterica"


## remove plasmids headers, concatenate with "n" break
#if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/fasta/",sp1,"/*.fna"))}
#if (!file.exists(paste0(dir,"data/external/fasta/",sp2,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/fasta/",sp2,"/*.fna"))}
#
## concatenate all files to one
#if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("cat ",dir,"data/external/fasta/",sp1,"/*.fna > ",dir,"data/external/fasta/",sp1,".fna"))}
#if (!file.exists(paste0(dir,"data/external/fasta/",sp2,".fna"))) {system(paste0("cat ",dir,"data/external/fasta/",sp2,"/*.fna > ",dir,"data/external/fasta/",sp2,".fna"))}
#
## find matches between two files
#soft <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/refseq_filter_download/bfmem/bfmem "
## soft <- "/home/misha/Documents/Development/Software/copmem2-main/copmem2 "

comm = paste0("mkdir -p ",dir,"data/processed/",sp2,"_",sp1,"/")
system(comm)

print("remove query headers with > (do only once and the comment out)")
comm2 =paste0("sed '/^>/d' ",dir,"data/processed/pairs/",sp2,"_",sp1,".tsv >",dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".tsv") 
system(comm2,wait=TRUE)

myFile = paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".tsv")

#if (file.size(myFile) > 0){
comp <- read.table(myFile, header = FALSE,  colClasses=c("character", "integer","integer","integer")); colnames(comp) <- c("header","start","foo","length")


comp = comp[which(comp[,4]>299),] # Need to add that otherwise h5 file is too big
comp <- comp[,c(1,2,4)]
#}else{
#comp = 	
#
#
#}


comp <- (comp %>% group_by_at(vars(header,start,length)) %>% summarise(weight = n(),, .groups="keep"))
write.table(comp,paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,"_weighted.tsv"),append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)

rTable <- wtd.table(as.integer(comp$length),w=as.numeric(comp$weight))
write.table(data.frame(r=as.integer(names(rTable)),m=as.integer(rTable)),paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".tsv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

comp <- data.table(start=as.integer(comp$start),end=as.integer(comp$start+comp$length-1),weight=as.numeric(comp$weight),header=comp$header)
comp <- split(comp[,1:3] , f = comp$header )
system(paste0("rm -f ",dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"),wait=TRUE)
h5save(comp, file=paste0(dir,"data/processed/",sp2,"_",sp1,"/",sp2,"_",sp1,".h5"))



###################################### 


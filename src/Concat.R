args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]

# dir <- "/home/misha/Documents/Development/ThreeFold/"
#dir <- "/scratch/ws1/msheinman-msheinman/ThreeFold/"
dir <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/"
# sp1 <- "Escherichia_coli"
# sp2 <- "Klebsiella_pneumoniae"
# sp2 <- "Salmonella_enterica"


# remove plasmids headers, concatenate with "n" break
if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/fasta/",sp1,"/*.fna"))}

# concatenate all files to one
if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("cat ",dir,"data/external/fasta/",sp1,"/*.fna > ",dir,"data/external/fasta/",sp1,".fna"))}



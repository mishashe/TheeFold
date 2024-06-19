args = commandArgs(trailingOnly=TRUE)
sp1 <- args[1]
sp2 <- args[2]

# dir <- "/home/misha/Documents/Development/ThreeFold/"
#dir <- "/scratch/ws1/msheinman-msheinman/ThreeFold/"
dir <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/"
# sp1 <- "Escherichia_coli"
# sp2 <- "Klebsiella_pneumoniae"
# sp2 <- "Salmonella_enterica"


# remove plasmids headers, concatenate with "n" break
#if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/fasta/",sp1,"/*.fna"))}
#if (!file.exists(paste0(dir,"data/external/fasta/",sp2,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/fasta/",sp2,"/*.fna"))}
#
# concatenate all files to one
#if (!file.exists(paste0(dir,"data/external/fasta/",sp1,".fna"))) {system(paste0("cat ",dir,"data/external/fasta/",sp1,"/*.fna > ",dir,"data/external/fasta/",sp1,".fna"))}
#if (!file.exists(paste0(dir,"data/external/fasta/",sp2,".fna"))) {system(paste0("cat ",dir,"data/external/fasta/",sp2,"/*.fna > ",dir,"data/external/fasta/",sp2,".fna"))}

# find matches between two files
soft <- "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/ThreeFold/src/bfmem/bfmem "
# soft <- "/home/misha/Documents/Development/Software/copmem2-main/copmem2 "
system(paste0("mkdir -p ",dir,"data/processed/",sp2,"_",sp1,"/"))
Command <- paste0(soft," -s b -t 1 -l 300  -o ",dir,"data/processed/pairs/",sp2,"_",sp1,".tsv -r ",dir,"data/external/fasta/",sp1,".fna"," -q ",dir,"data/external/fasta/",sp2,".fna")
system(Command)


dir <- "/home/misha/Documents/Development/ThreeFold/"
sp1 <- "Escherichia_coli"
sp2 <- "Klebsiella_pneumoniae"
# sp2 <- "Salmonella_enterica"


# remove plasmids headers, concatenate with "n" break
if (!file.exists(paste0(dir,"data/external/",sp1,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/",sp1,"/*.fna"))}
if (!file.exists(paste0(dir,"data/external/",sp2,".fna"))) {system(paste0("sed -i -r '1 ! s/^(>[^ ]+) .*/n/' ",dir,"data/external/",sp2,"/*.fna"))}

# concatenate all files to one
if (!file.exists(paste0(dir,"data/external/",sp1,".fna"))) {system(paste0("cat ",dir,"data/external/",sp1,"/*.fna > ",dir,"data/external/",sp1,".fna"))}
if (!file.exists(paste0(dir,"data/external/",sp2,".fna"))) {system(paste0("cat ",dir,"data/external/",sp2,"/*.fna > ",dir,"data/external/",sp2,".fna"))}

# find matches between two files
soft <- "/home/misha/Documents/Development/Software/bfmem/bfmem "
# soft <- "/home/misha/Documents/Development/Software/copmem2-main/copmem2 "
system(paste0("mkdir -p ",dir,"data/processed/",sp2,"_vs_",sp1,"/"))
Command <- paste0(soft," -s b -t 5 -l 300  -o ",dir,"data/processed/",sp2,"_vs_",sp1,"/",sp2,"_vs_",sp1,".tsv -r ",dir,"data/external/",sp1,".fna"," -q ",dir,"data/external/",sp2,".fna")
system(Command)




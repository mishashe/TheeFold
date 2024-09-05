



binSize = 0.15
# delta = 0.55 


coarse_grain <- function(m, rB,  binSize = binSize){
  
  if (binSize == 0){
    rMax = max(rB+1)
    if (rMax>5*10^4){ rMax = 5*10^4}
      rBnew = 1:rMax
  }else{
    # rBnew <- c(10^seq(log10(1)+0.1,5,binSize))
    rBnew <- c(seq(0,36,3),10^seq(log10(36)+0.1,5,binSize))
  }
  mnew <- rep(0,length(rBnew)-1)
  diffrB <- diff(rB)
  diffrBnew <- diff(rBnew)
  for (i in 2:length(rBnew))
  {
    Ind <- which(rB<=rBnew[i] & rB>rBnew[i-1])
    mnew[i-1] <- sum(m[Ind])/diff(rBnew)[i-1]
  }
  Mids = (rBnew[-1] + rBnew[-length(rBnew)])/2
  MLD = cbind( Mids, mnew)
  return(MLD)
}



allFiles = list.files("~/HGTnew/multi_comparisons/data/mobilome/alignments/pairs/")
allFiles = allFiles[grep (allFiles,pattern = "lastZ",invert =T)]
allFiles = allFiles[grep (allFiles,pattern = "mld",invert =T)]             
x = 1:10000

for (i in allFiles){
pdf(paste0("~/HGTnew/multi_comparisons/plots/mobilome/new_plots/",i,".pdf"))

  MLD_file = paste0("~/HGTnew/multi_comparisons/data/mobilome/alignments/pairs/",i,".mld")
  if (file.exists(MLD_file)){
    print("MLD exist")
  }else{
  comm = paste0("cut -f 4 ~/HGTnew/multi_comparisons/data/mobilome/alignments/pairs/",i," | sort -n |uniq -c |grep -v '>'")
  a = read.table(pipe(comm))
  write.table(a, file = MLD_file, append = FALSE, quote = FALSE, sep = "\t",row.names = F,col.names = F)
  }
  
  a = read.table(MLD_file,stringsAsFactors = FALSE)

  MLD2 = coarse_grain(m  = a$V1, rB = a$V2,binSize=0.15)

  plot(MLD2, xlim = c(1,max(a$V2)*5),log="xy", pch = "+", main = i)
myTwenty = which.min(abs(MLD2[,1]-20))

  lines(x,x^-3 * 20^3 * MLD2[myTwenty,2],lty = 2, col = 2)
  lines(x,x^-4 * 20^4 * MLD2[myTwenty,2],lty = 2, col = 4)

  print(i)
  print(date())
  dev.off()
}


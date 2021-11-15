####:: removeT7dups.R ::####
#############################

####################################################################
######:: Remove T7 duplicates -- based on read1 start position ::#######
####################################################################

setwd("~/sciTIP_pool/beds_rmdup/")

system("find . -type f -size 0 -delete") ## remove cells with 0 bytes to avoid error below   
allBEDs=list.files(pattern=".bed")
for (f in allBEDs){
  tip=read.delim(f, header=F)
  tip2=tip[,c(1,2,3)] ### will need to change this is using bedpe input --> tip2=tip[,c(1,2,6)]
  names(tip2)=c("CHR","StartR1","EndR2")
  tip2$ID1=paste(tip2$CHR,tip2$StartR1,sep="_")
  tip2$ID2=paste(tip2$V1,tip2$EndR2,sep="_")
  tip2$index=1:nrow(tip2)
  undups=tip2[!duplicated(tip2$ID1),c("index")]
  tip3=tip[undups,]
  write.table(tip3, file = paste0(f,"_rmT7dup.bed"), row.names = FALSE, quote = FALSE, sep = "\t", col.names=FALSE)
  cat(paste("Iteration", f, "was finished.\n"))
}

system("for i in *bed_rmT7dup.bed; do mv $i ${i%bed_rmT7dup.bed}rmT7dup.bed; done")
system("mkdir ../beds_rmT7dup")
system("mv *.rmT7dup.bed ../beds_rmT7dup/")

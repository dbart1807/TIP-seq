###:: deindex_sciTIP_Novaseq.R ::###
####################################
####################################
# Authored by Daniel Bartlett and Vishnu Dileep

options(mc.cores=20) ## <- set number of computing cores to use
setwd("~/sciTIP_pool/") ## <- set to directory containing pool.fastq files and barcode files.

####:: In preparation for deindexing ::#########

# 1. Prepare barcode.txt files (tab delimited)
#    If multiple samples/experimental treatments were multiplexed, edit the 'sample' column in r5_barcodes.csv file to reflect sample name - this will be appended to the output fastq file. 
#    Note that deindexing is slower the more possible index combinations there are to search,
#    so to save computation time user should delete any unused indexes from barcode.txt files.

# Note: Illumina sequencing platforms that use the reverse complement indexing workflow (ie. Novaseq v1.5 reagents) 
#       will need the reverse complement r5 and i5 index files (included on GitHub as: i5_barcodes.revComp.txt; r5_barcodes.revComp.txt)
#       and will also flip the order in which the indexes are found in the Read1 index sequence header.
#       *If using Illumina forward strand index chemistry, set script to pull from different nucleotides in read 1 index file
#       (more info at https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-08.pdf)

# 2. check that I2 index is 29 bp long and is appended to R1/R2 fastq read headers
system("head -c 2000 Undetermined_S0_L001_R1_001.fastq.gz | zcat 2>/dev/null | head -10")    # head gzipped file.

# 3. concatenate Lane1 and Lane2 files 
system("for i in *L001*; do cat $i ${i/L001/L002} > ${i/_L001/}; done")

# 4. rename files # note: cat command above unzipped files so remove .gz from names
system("mv Undetermined_S0_I1_001.fastq.gz sciTIP_pool_I1_001.fastq")
system("mv Undetermined_S0_R1_001.fastq.gz sciTIP_pool_R1_001.fastq")
system("mv Undetermined_S0_R2_001.fastq.gz sciTIP_pool_I2_001.fastq") ## Index2 = *_R2.fastq file in this novaseq bcl2fastq scheme
system("mv Undetermined_S0_R3_001.fastq.gz sciTIP_pool_R2_001.fastq") ## Read2 = *_R3.fastq file in this novaseq bcl2fastq scheme

# 5. split up fastq files into 50M read chunks (better memory usage)
system("split -l 200000000 sciTIP_pool_R1_001.fastq sciTIP_pool_R1_001.fastq.split.") # splits fastq into 50M read chunks (to limit memory usage)
system("split -l 200000000 sciTIP_pool_R2_001.fastq sciTIP_pool_R2_001.fastq.split.") # splits fastq into 50M read chunks (to limit memory usage)

library(data.table)
library(Biostrings)
library(parallel)
library(tidyverse)

################################
###:: read barcode files ::#####
################################

i7_index=read.delim("i7_barcodes.txt", header=T)
i5_index=read.delim("i5_barcodes.revComp.txt", header=T)
r5_index=read.delim("r5_barcodes.revComp.txt", header=T)

####:: make all possible indices with format i5+i7+r5 ::#######
i7_index$i7_comb=paste(i7_index$i7_sequence,i7_index$i7_name, sep=",")
i5_index$i5_comb=paste(i5_index$i5_sequence,i5_index$i5_name, sep=",")
r5_index$r5_name=paste(r5_index$sample,r5_index$r5_name, sep="_")
r5_index$r5_comb=paste(r5_index$r5_sequence,r5_index$r5_name, sep=",")
all_combs <- expand.grid(i5_index$i5_comb, i7_index$i7_comb, r5_index$r5_comb)

####:: separate index_names from index_sequence ::##
my_split_fun1=function(x)  unlist(strsplit(as.character(x),split=","))[c(F,T)]
my_split_fun2=function(x)  unlist(strsplit(as.character(x),split=","))[c(T,F)]
all_combs=cbind(apply(all_combs,2,my_split_fun1),
                apply(all_combs,2,my_split_fun2))
all_combs=as.data.frame(all_combs)

#####:: combine all indexes combinations
all_combs$names=paste(all_combs[,3],all_combs[,1],all_combs[,2],sep="_")
all_combs$combindex=paste0(all_combs[,4],all_combs[,5],all_combs[,6])
all_combs=as.data.frame(all_combs)
all_combs[,1:6]=NULL
all_combs=as.data.frame(all_combs)

###################################################################
####:: Iterate deindexing over all split R1/R2 fastq files ::######
###################################################################

chunk = unique(unlist(strsplit(list.files(pattern = ".fastq.split"), '_001.fastq.'))[c(FALSE,TRUE)]) # number of iterations (ie. split fastqs)
length(chunk) # number of iterations (ie. split fastq files)

system("mkdir deindexed_fastq")
for(i in 1:length(chunk)){
  
  ###:: load in R1/R2 fastq files ::#####
  F1=fread(paste0("sciTIP_pool_R1_001.fastq.",chunk[i]), header=FALSE,col.names = "L",sep="") 
  F2=fread(paste0("sciTIP_pool_R2_001.fastq.",chunk[i]), header=FALSE,col.names = "L",sep="") 
  
  ###:: Extract i5 and r5 indexes from Index2 and recombine all indices i5+i7+r5  ::###  NOTE: If NOT using reverse complement chemistry, you will need to change these numbers... order will be reversed for i5/r5 index sequence.
  l=nrow(F1)
  rawIndex=F1$L[seq(1,l,4)]
  rawIndex=strsplit(rawIndex,split=":")
  i7s=unlist(lapply(rawIndex,'[[',11))
  r5i5=unlist(lapply(rawIndex,'[[',8))
  i5s=unlist(lapply(r5i5,substring,22,29))
  r5s=unlist(lapply(r5i5,substring,1,6))
  
  #rm(rawIndex,r5i5)
  index_file=paste0(i5s,i7s,r5s)
  #rm(i5s,r5s,i7s)
  
  ####frequency of all index###
  t=table(index_file)
  t=stack(t)
  t=t[order(t$values,decreasing = T),]
  names(t)=c("Hits","combindex")
  head(t)
  
  ###:: true index combinations ::###
  all_truecombs=merge(all_combs,t,by="combindex")
  sum(all_truecombs$Hits)/length(index_file) ###percentage of total reads that have true index##
  all_truecombs=subset(all_truecombs, Hits>=50) ###only indexes with more than X frequency##
  all_truecombs=all_truecombs[order(all_truecombs$Hits,decreasing = T),] ### sorts by hits
  
  ###################################
  ####:: deindexing functions ::#####
  
  Ori.L.read=function(x) 4*x-2
  
  splitbyBC=function(f) {
    
    pmatch= index_file == all_truecombs$combindex[f]
    pmatch.count=sum(pmatch)
    pmatch.rows=which(pmatch)
    
    ####:: convert to original fastq rownumbers ::####
    pmatch.ori.rows=Ori.L.read(pmatch.rows)
    fullrows=sort(sapply(pmatch.ori.rows,function(x) seq(x-1,x+2,1)))
    rm(pmatch,pmatch.rows,pmatch.ori.rows)
    
    ####:: write read 1 :####
    match.fastq1=F1[fullrows,"L"]
    matched.filename1=paste0("./deindexed_fastq/",all_truecombs$names[f],"_R1",".fastq.",chunk[i])
    fwrite(match.fastq1,matched.filename1,sep="\t", row.names = FALSE,col.names = FALSE, quote = FALSE)
    
    ####:: write read 2 ::####
    match.fastq2=F2[fullrows,"L"]
    matched.filename2=paste0("./deindexed_fastq/",all_truecombs$names[f],"_R2",".fastq.",chunk[i])
    fwrite(match.fastq2,matched.filename2,sep="\t",row.names=FALSE,col.names = FALSE, quote = FALSE)
    
    rm(fullrows,match.fastq1,match.fastq2)
    return(pmatch.count)
  }
  
  ##finding matches and writing parsed fastqs###
  out=unlist(mclapply(1:nrow(all_truecombs),splitbyBC))
  head(out)
  write.table(all_truecombs,file = paste0("All_index_combinations.",chunk[i],".txt"),sep="\t",row.names=FALSE,quote = FALSE)
  
  cat(paste("Fastq chunk ", i, "of",length(chunk),"was finished.\n"))
}

####:: combine split All_index_combinations.txt files ::####
statfiles = lapply(list.files(pattern = "All_index_combinations.split"), function(x)read.table(x, header=T)) 
stat=statfiles %>% reduce(full_join, by = c("names","combindex"))
stat$Hits_sum=rowSums(stat[,3:ncol(stat)])
comb_stat=stat[,c('names','combindex','Hits_sum')]
write.table(comb_stat,file = "stats_deindexing.txt",sep="\t",row.names=FALSE,quote = FALSE)
system("rm All_index_combinations*")

####:: recombine all fastq chunks ::####
setwd(paste0(getwd(),"/deindexed_fastq/"))
fqcomb = unique(unlist(strsplit(list.files(pattern = ".fastq.split."), '.fastq.'))[c(TRUE,FALSE)]) # number of iterations (ie. split fastqs)
length(fqcomb) # number of iterations (ie. split fastq files)
for(i in 1:length(fqcomb)){
  cmd<-paste0("cat ",fqcomb[i],"* > ",fqcomb[i],".fastq")
  system(cmd)
}


#####################
####:: cleanup ::#### 
### enable cleanup below when script completed successfully.


#system("gzip ../*fastq") # gzip original pooled fastq files
#system("rm ../*_001.fastq.split.*") # delete pooled split.fastq files
#system("rm *.fastq.split.*") # delete deindexed split.fastq files

#### Next: run "mapping.sciTIP.R" script

###:: mapping_mm10.sciTIP.R ::#####
####################################
### run after deindexing script
####################################

library(travis)
set.cores=20

btindex.path="/home/share/references/assembly/mm10"
chromsize.path="/home/share/references/chromsizes/mm10.chrom.sizes"

setwd("~/sciTIP_pool/deindexed_fastq/")

system("mkdir ../mapped_mm10")
system("mkdir ../stats_mapping")
fq <- unique(strsplit(list.files(pattern = "R1.fastq"), '_R1.fastq'))
length(fq) # number of cells 
for(i in 1:length(fq)){
  Sam=bowtie2(read1files = paste0(fq[i],"_R1.fastq"), read2files = paste0(fq[i],"_R2.fastq"), indexfile = btindex.path, 
            alignMode = "very-sensitive-local", appendIndexToName=TRUE, unaligned =FALSE, 
            reorder=TRUE, threads=set.cores, minInsertSize = 10, maxInsertSize = 700)
  cmd=paste0("mv ",fq[i],"_R1.fastq" )
  system("for i in *R1_mm10.sam; do mv $i ${i%_R1_mm10.sam}.mm10.sam; done")
  system("for i in *_R1_mm10.sam.log; do mv $i ${i%_R1_mm10.sam.log}.mm10.sam.log; done")
  system("mv *.mm10.sam ../mapped_mm10")
  system("mv *sam.log ../stats_mapping")
  }

print("bowtie2 mapping complete!")

#####:: cleanup ::#####
system("gzip *.fastq")

# Note: function above yields the following bowtie2 parameters (check sam header)
#    bowtie2-align-s --wrapper basic-0 --very-sensitive-local -p 28 --no-mixed --no-discordant 
#           --reorder -I 10 -X 700 -q -x /home/share/references/assembly/hg38 --passthrough 
#           -1 WT-k9m3_S508_N704_r5025_R1.fastq -2 WT-k9m3_S508_N704_r5025_R2.fastq


## Next: run "removeDups.R"

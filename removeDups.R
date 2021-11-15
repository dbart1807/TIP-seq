######:: removeDups.R ::##

### Run after mapping script

####################################################################
######:: Remove duplicates -- samtools fixmate and markdup ::#######
####################################################################

setwd("~/sciTIP_pool/mapped_mm10/")

##:: filter out <q10 reads and convert to name sorted BAM
nsortBam=paste0("for i in *mm10.sam; do samtools view -q 10 -@ ",set.cores," -u $i | samtools sort -n -o ${i%sam}q10.nsort.bam; done")
head (nsortBam)  ## delete
system(nsortBam)

### enable only if bam conversion successful!!
#system("rm *.sam")

##:: run fixmate
fixBam=paste0("for i in *q10.nsort.bam; do samtools fixmate -@ ",set.cores," -m $i ${i%nsort.bam}fixmate.bam; done")
system(fixBam)

##:: sort by coordinate
psortBam=paste0("for i in *fixmate.bam; do samtools sort -@ ",set.cores," -o ${i%bam}psort.bam $i; done")
system(psortBam)

##:: run markdup
mkdupBams=paste0("for i in *fixmate.psort.bam; do samtools markdup -r -@ ",set.cores," $i ${i%mate.psort.bam}mkdup.bam; done")
system(mkdupBams)

# name sort (for bamToBed conversion)
nsort_mkdupBams=paste0("for i in *.q10.fixmkdup.bam; do samtools sort -n -@ ",set.cores," -o ${i%bam}nsort.bam $i; done")
system(nsort_mkdupBams)

print("samtools fixmate & markdup complete!")

##:: bamToBed conversion
nsort_mkdupBams= list.files(pattern = "*mkdup.nsort.bam")
Beds=bamToBed(nsort_mkdupBams, paired = TRUE, threads=set.cores, sortThreads = set.cores)  # requires nsort bam

print("BAM to BED conversion complete!")

system("mkdir ../beds_rmdup")
system("for i in *nsort.bed; do mv $i ../beds_rmdup/${i%nsort.bed}bed; done")

### CLEANUP: enable only if fixmate.mkdup conversion successful!! to clear disk space if necessary
system("rm *.q10.fixmate*")
system("rm *.q10.fixmkdup.bam")


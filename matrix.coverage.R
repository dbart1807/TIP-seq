####:: matrix.coverage.R ::####

##################################################
####:: calculate coverage and build matrix ::#####
##################################################

library(travis)
set.cores=20
build = "mm10"
bin_size=5000
step_size=5000
scalar=1

chromsize_path.mm10="/home/share/references/chromsizes/mm10.chrom.sizes"

setwd("~/sciTIP_pool/beds_rmdup/")  # path to beds
allBEDs <- list.files(pattern="fixmkdup.bed")

windows <- bedtoolsMakeWindows(bedfiles = chromsize_path.mm10, windowsize=bin_size,stepsize=step_size, threads = set.cores, genome=TRUE,
                               outnames = paste0("~/bin/",build,"_win",bin_size,".step",step_size,".bed"))

system("for i in *bed; do bedSort $i $i; done")
cov=paste0("for i in *.bed; do bedtools coverage -counts -a ",windows," -b $i > ${i%bed}win5k.bg; done")
system(cov)

###:: bigwig conversion (optional, for IGV)  # may need to run in command line...
bg2bw=paste0("for j in *.bg; do bedGraphToBigWig $j ",chromsize_path.mm10," ${j%bg}bw; done")
system(bg2bw) 

system("mkdir ../coverage_win5k")
system("mv *.win5k.bg ../coverage_win5k")

system("mkdir ../coverage_win5k/win5k_bigwigs")
system("mv *.win5k.bw ../coverage_win5k/win5k_bigwigs")

###########################################################
########--- build single-cell coverage matrix --- #########
###########################################################
set.cores= 20
setwd("~/sciTIP_pool/coverage_win5k/") 
allBgs <- list.files(pattern="win5k.bg")

a= bgRead(allBgs,threads=set.cores)
coords= read.tsv(allBgs[1])
fdata=cbind(coords[,1:3],a)
nam=names(fdata)
nam=nam[4:length(nam)]
nam=unlist(strsplit(nam,".mm10.q10"))[c(TRUE,FALSE)]
nam=c("CHR","START","END",nam)
names(fdata)=nam
write.table(fdata, file=gzfile("sciTIP_counts_matrix.txt.gz"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

##:: cleanup ::##
system("rm *.win5k.bg") # optional: delete bedgraphs, do not need the after creating matrix and making bigwigs.

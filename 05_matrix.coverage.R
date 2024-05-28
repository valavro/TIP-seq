####:: matrix.coverage.R ::####

##################################################
####:: calculate coverage and build matrix ::#####
##################################################

library(travis)
set.cores=8
build = "mm10"
bin_size=5000
step_size=5000
scalar=1

chromsize_path.mm10="/home/cjm933/datadir/TIP_scripts/Data/chromInfo.main.tsv"

setwd("/maps/projects/dan1/people/cjm933/TIP_scripts/Data/beds_rmT7dup")  # path to beds
allBEDs <- list.files(pattern=".bed")

windows <- bedtoolsMakeWindows(bedfiles = chromsize_path.mm10, windowsize=bin_size,stepsize=step_size, threads = set.cores, genome=TRUE,
                               outnames = paste0("/home/cjm933/datadir/TIP_scripts/Data/",build,"_win",bin_size,".step",step_size,".bed"))

#system("for i in *bed; do bedSort $i $i; done")

# cov=paste0("for i in *sorted.bed; do bedtools coverage -counts -a ",windows," -b $i > ${i%bed}win5k.bg; done")
# system(cov)

for(file in allBEDs){
  bedSort(file, threads=set.cores)
  cat(paste0("Done with file ", file))
}

cov=paste0("for i in *.bed; do bedtools coverage -counts -a ",windows," -b $i > ${i%bed}win5k.bedgraph; done")
system(cov)


allBedGraphs =list.files(pattern=".bedgraph")


###:: bigwig conversion (optional, for IGV)  # may need to run in command line...
#bg2bw=paste0("for j in *.bg; do bedGraphToBigWig $j ",chromsize_path.mm10," ${j%bg}bw; done")
#system(bg2bw) 

for(file in allBedGraphs){
  bigwigfile <- bedGraphToBigWig(file, chromsizes = chromsize_path.mm10, threads=set.cores)
  save(bigwigfile, file = paste0(file, ".bigwig"))
  print(paste0("Done with file ", file))
}

system("for i in *bedgraph.bigwig; do mv $i ${i%bedgraph.bigwig}bigwig; done")



system("mkdir ../coverage_win5k")
system("mv *.win5k.bedgraph ../coverage_win5k")

system("mkdir ../coverage_win5k/win5k_bigwigs")
system("mv *.bigwig ../coverage_win5k/win5k_bigwigs")

###########################################################
########--- build single-cell coverage matrix --- #########
###########################################################
set.cores= 8
setwd("/maps/projects/dan1/people/cjm933/TIP_scripts/Data/coverage_win5k/") 
allBgs <- list.files(pattern="win5k.bedgraph")

a= bgRead(allBgs[7:11],threads=set.cores)
coords= read.tsv(allBgs[1])
fdata=cbind(coords[,1:3],a)
nam=names(fdata)
nam=nam[4:length(nam)]
nam=unlist(strsplit(nam,".genome.q10"))[c(TRUE,FALSE)]
nam=c("CHR","START","END",nam)
names(fdata)=nam
write.table(fdata, file=gzfile("H3K27me3_counts_matrix.txt.gz"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

##:: cleanup ::##
system("rm *.win5k.bg") # optional: delete bedgraphs, do not need the after creating matrix and making bigwigs.

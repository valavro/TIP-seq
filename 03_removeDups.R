######:: removeDups.R ::##

### Run after mapping script

####################################################################
######:: Remove duplicates -- samtools fixmate and markdup ::#######
####################################################################

setwd("/maps/projects/dan1/people/cjm933/TIP_scripts/Data/mapped_mm10")
set.cores=8


##:: filter out <q10 reads and convert to name sorted BAM
nsortBam=paste0("for i in *genome.sam; do samtools view -q 10 -@ ",set.cores," -u $i | samtools sort -n -o ${i%sam}q10.nsort.bam; done")
head (nsortBam)  ## delete
#system("module load samtools/1.15.1")
system(nsortBam)
#samtools view sampleA_r5001_S502_N712.sam -q 10 -@ 8 -u | samtools sort -n -o sampleA_r5001_S502_N712.q10.nsort.bam





### enable only if bam conversion successful!!
#system("rm *.sam")




##:: run fixmate

fixBam=paste0("for i in *q10.nsort.bam; do samtools fixmate -@ ",set.cores," -m $i ${i%nsort.bam}fixmate.bam; done")
system(fixBam)
#samtools fixmate sampleA_r5001_S502_N712.q10.nsort.bam -@ 8 -m sampleA_r5001_S502_N712.q10.fixmate.bam




##:: sort by coordinate

psortBam=paste0("for i in *fixmate.bam; do samtools sort -@ ",set.cores," -o ${i%bam}psort.bam $i; done")
system(psortBam)
#samtools sort sampleA_r5001_S502_N712.q10.fixmate.bam -@ 8 -o sampleA_r5001_S502_N712.q10.fixmate.psort.bam





##Option1:: Run markdup to deduplicate based on both reads

mkdupBams=paste0("for i in *fixmate.psort.bam; do samtools markdup -r -@ ",set.cores," $i ${i%mate.psort.bam}mkdup.bam; done")
system(mkdupBams)
#samtools markdup sampleA_r5001_S502_N712.q10.fixmate.psort.bam -r -@ 8 sampleA_r5001_S502_N712.q10.fixmate.psort.mkdup.bam




# ### Option2:: 
# ##Separate R1s 
# 
# getR1s=paste0("for i in *fixmate.psort.bam; do samtools view -b -f 64 -@ ",set.cores," $i -o ${i%mate.psort.bam}mate.psortR1.bam; done")
# system(getR1s)
# 
# ## Remove duplicates and 
# 
# mkdupBamsR1=paste0("for i in *fixmate.psortR1.bam; do samtools markdup -r -@ ",set.cores," $i ${i%mate.psortR1.bam}mkdupR1.bam; done")
# system(mkdupBamsR1)
# 
# ## Filter paired BAM 
# 
# # Step 1: Convert R1-only BAM to list of read names
# system("samtools view -h *N712*dupR1.bam | awk '$1 ~ /^@/ {print; next} {print $1}' > read_namesN712.txt")
# 
# # Step 2: Filter paired read BAM using the list of read names
# system("samtools view -h sampleA_r5001_S502_N712.genome.q10.fixmkdup.nsort.bam | grep -Fwf read_namesN712.txt - | samtools view -b - > filtered_pairedN712.bam")
# system("samtools sort -@ 8 filtered_pairedN712.bam -o filtered_pairedN712.psort.bam")
# 
# # Index the filtered BAM
# system("samtools index filtered_pairedN712.psort.bam")
# system("samtools sort sampleA_r5001_S502_N712.genome.q10.fixmkdup.nsort.bam -o N712.psort.bam")
# system("samtools sort sampleA_r5001_S502_N712.genome.q10.fixmkdupR1.bam -o N712R1.psort.bam")
# system("samtools index N712.psort.bam")








# name sort (for bamToBed conversion)
nsort_mkdupBams=paste0("for i in *.q10.fixmkdup.bam; do samtools sort -n -@ ",set.cores," -o ${i%bam}nsort.bam $i; done")
system(nsort_mkdupBams)
#samtools sort sampleA_r5001_S502_N712.q10.fixmate.psort.mkdup.bam -n -@ 8 -o sampleA_r5001_S502_N712.q10.fixmkdup.nsort.bam









print("samtools fixmate & markdup complete!")

##:: bamToBed conversion
nsort_mkdupBams= list.files(pattern = "*mkdup.nsort.bam")
Beds=bamToBed(nsort_mkdupBams, paired = TRUE, threads=set.cores, sortThreads = set.cores)  # requires nsort bam
#bedtools bamtobed -i sampleA_r5001_S502_N712.q10.fixmkdup.nsort.bam > sampleA_r5001_S502_N712.bed


print("BAM to BED conversion complete!")

system("mkdir ../beds_rmdup")
system("for i in *nsort.bed; do mv $i ../beds_rmdup/${i%nsort.bed}bed; done")

### CLEANUP: enable only if fixmate.mkdup conversion successful!! to clear disk space if necessary
system("rm *.q10.fixmate*")
system("rm *.q10.fixmkdup.bam")


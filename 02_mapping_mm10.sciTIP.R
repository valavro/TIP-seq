###:: mapping_mm10.sciTIP.R ::#####
####################################
### run after deindexing script
####################################

library(travis)
set.cores=8

btindex.path="/maps/projects/dan1/data/GENOMICS/DATABASES/GENOME/Mouse/mm10/Bowtie2Index/genome"
chromsize.path="/home/cjm933/datadir/TIP_scripts/Data/chromInfo.main.tsv"

setwd("/maps/projects/dan1/people/cjm933/TIP_scripts/Data/deindexed_fastq")

system("mkdir ../mapped_mm10")
system("mkdir ../stats_mapping")
fq <- unique(strsplit(list.files(pattern = "R1.fastq"), '_R1.fastq'))
length(fq) # number of cells 
for(i in 7:length(fq)){
  cat(i)
  Sam=bowtie2(read1files = paste0(fq[i],"_R1.fastq"), read2files = paste0(fq[i],"_R2.fastq"), indexfile = btindex.path, 
            alignMode = "very-sensitive-local", appendIndexToName=TRUE, unaligned =FALSE, 
            reorder=TRUE, threads=set.cores, minInsertSize = 10, maxInsertSize = 700)
  cmd=paste0("mv ",fq[i],"_R1.fastq" )
  system("for i in *R1_genome.sam; do mv $i ${i%_R1_genome.sam}.genome.sam; done")
  system("for i in *_R1_genome.sam.log; do mv $i ${i%_R1_genome.sam.log}.genome.sam.log; done")
  system("mv *.genome.sam ../mapped_mm10")
  system("mv *sam.log ../stats_mapping")
  }

print("bowtie2 mapping complete!")

#####:: cleanup ::#####
system("gzip *.fastq")

# Note: function above yields the following bowtie2 parameters (check sam header)
#    bowtie2-align-s --wrapper basic-0 --very-sensitive-local -p 28 --no-mixed --no-discordant 
#           --reorder -I 10 -X 700 -q -x /home/share/references/assembly/hg38 --passthrough 
#           -1 WT-k9m3_S508_N704_r5025_R1.fastq -2 WT-k9m3_S508_N704_r5025_R2.fastq

### I have removed the --passthrough flag because it led to "@" at the beginning of lines and I have redirected the outputs and logs
#bowtie2-align-s --wrapper basic-0 --very-sensitive-local -p 28 --no-mixed --no-discordant --reorder -I 10 -X 700 -q -x /maps/projects/dan1/data/GENOMICS/DATABASES/GENOME/Mouse/mm10/Bowtie2Index/genome -1 sampleA_r5001_S502_N712_R1.fastq -2 sampleA_r5001_S502_N712_R2.fastq -S ../mapped_mm10/sampleA_r5001_S502_N712.sam > ../mapped_mm10/bowtie_log.txt 2>&1


## Next: run "removeDups.R"

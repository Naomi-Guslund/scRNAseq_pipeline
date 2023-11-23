#From_fastq_to_Digital_expression_matrix

#Make sure all genome files are present. Descriptions of files needed are found in ames Nemesh, McCarroll Lab
#Drop-seq core computational protocol
#V2.0.0; September 28 2018


# Example fastq files from sequencing centre
## Sample_1_R1_001.fastq.gz
## Sample_1_R2_001.fastq.gz          

#Unzip files
gunzip *fastq.gz

#Run a slurm script


##################################### 
#!/bin/bash
##This is the master script for submitting jobs to preprocess reads to the grid. The time below is 48 hours
##change the name below
#SBATCH --job-name=1-COD-Sample_1
#SBATCH --account=nn9244k
#SBATCH --time=8:0:0
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=30G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#set -o errexit  # Recommended for easier debugging
## Load your modules
#module purge   # Recommended for reproducibility
module --force purge
module load StdEnv
module load FastQC/0.11.8-Java-1.8
module load Drop-seq_tools/2.3.0
module load Java/1.8.0_212
module load picard/2.18.27-Java-1.8
module load STAR/2.7.3a-GCC-8.3.0

##convert fastq files to bam files using picard-tools

java -jar $EBROOTPICARD/picard.jar \
FastqToSam \
F1=Sample_1_R1_001.fastq \
F2=Sample_1_R2_001.fastq \
O=Sample_1_raw_sorted.bam \
SORT_ORDER=queryname \
SAMPLE_NAME=Sample_1 \
2>Sample_1_fastq_to_bam_sorted_queryname.err

##1.a. Tag cell barcodes
TagBamWithReadSequenceExtended \
SUMMARY=Sample_1_tagcellbarcodes.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=false \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=Sample_1_raw_sorted.bam \
OUTPUT=Sample_1_raw_sorted.tagcellbarcode.bam \
2>Sample_1_tagcellbarcode.err

##1.b. Tag molecular barcodes
TagBamWithReadSequenceExtended \
SUMMARY=Sample_1_tagmolbarcodes.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=true \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=Sample_1_raw_sorted.tagcellbarcode.bam \
OUTPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode.bam \
2>Sample_1_tagmolbarcode.err


FilterBam \
TAG_REJECT=XQ \
INPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode.bam \
OUTPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt.bam \
2>Sample_1_filterbam.err


##1.c. Trim 5’ primer sequence
TrimStartingSequence \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5 \
INPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt.bam \
OUTPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim.bam \
OUTPUT_SUMMARY=Sample_1_5trim_report.txt \
2>Sample_1_trim5end.err


##1.d. Trim 3’ polyA sequence
PolyATrimmer \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true \
INPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim.bam \
OUTPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim_Atrim.bam \
OUTPUT_SUMMARY=Sample_1_Atrim_report.txt \
2>Sample_1_trimAtail.err


##1.e. convert SAM ¬file to Fastq file
java -Xmx4g -jar $EBROOTPICARD/picard.jar SamToFastq \
INPUT=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim_Atrim.bam \
FASTQ=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim_Atrim.fastq \
2>Sample_1_convertAtrimtofastq.err 


######1.f. alignment (index should be created in genomedir before running this step)
##mapping the reads from the bead tag + umi tag sorted + adapter trimmed + polyA trimmed files to the genome index for counting
##STAR --runMode genomeGenerate \
##--genomeFastaFiles GCF_902167405.1_gadMor3.0_nu_mt_for_Naomi_clean.fasta \
##--runThreadN 10 \
##--genomeSAindexNbases 13 \
##--sjdbGTFfile GCF_902167405.1_gadMor3.0_genomic.renamed.dropseq.final.gata2.gtf \
##--sjdbGTFtagExonParentTranscript Parent \
##--sjdbOverhang 74 \
##--genomeDir .

STAR \
--genomeDir . \
--runThreadN 10 \
--outFileNamePrefix star_Sample_1_ \
--readFilesIn Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim_Atrim.fastq \
2>Sample_1_star.err


##1.g. Sort STAR alignment in queryname order (STAR does not necessarily emit reads in the same order as the input)

java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=1 -Xmx4g -jar $EBROOTPICARD/picard.jar \
SortSam \
INPUT=star_Sample_1_Aligned.out.sam \
OUTPUT=star_Sample_1_Aligned.out.sorted.bam \
SORT_ORDER=queryname \
TMP_DIR=. \
2>Sample_1_star_filesort.err


##1.h. Merge STAR alignment tagged SAM to recover cell/molecular barcodes(STAR does not keep the tags in its result)
# make sure dictionary is present, made with picard. 



java  -Xmx4g -XX:ParallelGCThreads=10 \
-jar $EBROOTPICARD/picard.jar \
MergeBamAlignment \
REFERENCE_SEQUENCE=GCF_902167405.1_gadMor3.0_nu_mt_for_Naomi_clean.fasta \
UNMAPPED_BAM=Sample_1_raw_sorted.tagcellbarcode.tagmolbarcode_filt_5trim_Atrim.bam \
ALIGNED_BAM=star_Sample_1_Aligned.out.sorted.bam \
OUTPUT=star_Sample_1_Aligned.out.sorted.retagged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false \
2>Sample_1_mergealignment.err

##1.i TagReadWithGeneExon 

TagReadWithGeneFunction \
I=star_Sample_1_Aligned.out.sorted.retagged.bam \
O=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged.bam \
ANNOTATIONS_FILE=GCF_902167405.1_gadMor3.0_genomic.renamed.dropseq.final.gata3.gtf \
1> star_Sample_1-PreImm_exontagging_equal.out 2>star_Sample_1-PreImm_exontagging_equal.err


##Number of predicted cells =  700 . Decided to include 3x the number of barcodes to give large margin of error for including maximum information
##1.j. DetectBeadSynthesisErrors
DetectBeadSynthesisErrors \
I=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged.bam \
O=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean.bam \
OUTPUT_STATS=Sample_1-PreImm_cleaning_stats.txt \
SUMMARY=Sample_1-PreImm_cleaning_summary.txt \
NUM_BARCODES=2100 \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
2>Sample_1-PreImm_cleaning.err


###End of Alignment
##At this point, the alignment is completed, and your raw reads have been changed from paired reads
##to single end reads with the cell and molecular barcodes extracted, cleaned up, aligned, and prepared
##for DGE extraction.
##2.cell selection
##2.a create read counts for Selecting cells 

BamTagHistogram \
I=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean.bam \
O=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean_cell_readcounts.txt.gz \
TAG=XC

####################################################################
sbatch run_COD-Sample_1.sh

## create a cumulative distribution plot in R
module load R/3.6.2-fosscuda-2019b
R
options(scipen=2000)

a1=read.table("star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean_cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
x=cumsum(a1$V1)
x=x/max(x)
pdf("Sample_1-PreImm_bamtaghisto.pdf")
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
ylab="cumulative fraction of reads", xlim=c(1,160000))
abline(v=10000,col="black")
dev.off()
barcodes_Sample_1 <- head(a1, n=10000)
write.table(barcodes_Sample_1, file="barcodes_Sample_1.txt", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(barcodes_Sample_1[,2], file="barcodes_Sample_1_clean.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

##quit R
q()

## Create a digital expression matrix ready to be loaded into RStudio
module load Drop-seq_tools/2.3.0

DigitalExpression \
I=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean.bam \
O=star_Sample_1-PreImm_Aligned.out.sorted.retagged_gene_exon_tagged_clean.bam.dge.txt.gz \
SUMMARY=Sample_1-PreImm_gene_exon_tagged.dge.summary.txt \
CELL_BC_FILE=barcodes_Sample_1_clean.txt \
TMP_DIR=/cluster/work/users/naomiac/tmp


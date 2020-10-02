#!/bin/bash
#SBATCH -A b1025
#SBATCH -p buyin
#SBATCH -t 24:00:00
#SBATCH -m a
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1025/sdi0596/whole_genome_Exome_practice/scripts/
#SBATCH -o "%x.o%j"
#SBATCH --job-name=ERR2210079
#SBATCH --nodes=1
#SBATCH -n 15

##WES Tumor
#If there are any modules loaded, remove them.
module purge

export PATH=$PATH:/projects/p20742//tools/bin/
module load samtools/1.6
module load boost/1.56.0
module load gcc/4.8.3
module load bwa/0.7.12
module load vcftools/0.1.17
module load bamtools
# Setting up trimming for single end reads.




# Trim PE poor quality sequence with TRAILING:30 MINLEN:20 (see Trimmomatic documentation)
## phred64 all reads are getting dropped
java -Xmx8G -jar /projects/p20742//tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE \
-threads 15 -phred33 \
-trimlog /projects/b1025/sdi0596/whole_genome_Exome_practice/fastq/ERR2210079.trim.log \
/projects/b1025/sdi0596/whole_genome_Exome_practice/fastq/ERR2210079_R1.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/fastq/ERR2210079_R2.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R1.passed.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R1.not_passed.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R2.passed.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R2.not_passed.fastq.gz \
LEADING:10 TRAILING:10 SLIDINGWINDOW:5:0 MINLEN:25


### TRIMMED READS LOOK OKAY

bwa index -a bwtsw /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68

#### ALIGNMENT USING bwa. Using mm10 since dbSNP142 is buolt using mm10
bwa mem -t 15 /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-Y -K 10000000 -v 1 \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R1.passed.fastq.gz \
/projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079_R2.passed.fastq.gz \
 > /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.sam


 ##POST PROCESSING OF ALIGNED READS
 java -Xmx8G -jar /projects/b1025/sdi0596/tools/picard.jar CleanSam \
 --INPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.sam \
 --OUTPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.bam \
 -VALIDATION_STRINGENCY LENIENT

 samtools sort /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.bam -o /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.bam ## Worked






samtools faidx /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa


### wrong dbsnp being used plus no addgroup issue was the issue earlier
java -Xmx8G -jar /projects/b1025/sdi0596/tools/picard.jar AddOrReplaceReadGroups \
--INPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.bam \
--OUTPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.readgroups.bam \
--RGID 1 \
--RGLB Lib1 \
--RGPL ILLUMINA \
--RGPU Run1 \
--RGSM ERR2210079

java -Xmx8G -jar /projects/b1025/sdi0596/tools/picard.jar MarkDuplicates \
--INPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.readgroups.bam \
--OUTPUT /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.readgroups.marked.bam \
--METRICS_FILE /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.duplicate_metrics.txt  \
--REMOVE_DUPLICATES false --ASSUME_SORTED true \
-VALIDATION_STRINGENCY LENIENT

/projects/b1025/sdi0596/tools/gatk-4.1.3.0/gatk CreateSequenceDictionary -R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa


#### BASE RECALIBRATION###
/projects/b1025/sdi0596/tools/gatk-4.1.3.0/gatk BaseRecalibrator \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.readgroups.marked.bam \
--known-sites /projects/b1025/sdi0596/whole_genome_Exome_practice/dbsnp_data/mgp.v3.snps.rsIDdbSNPv137.vcf.gz \
--use-original-qualities \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.GATK4.pre.recal.table \


/projects/b1025/sdi0596/tools/gatk-4.1.3.0/gatk ApplyBQSR \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.cleaned.sorted.readgroups.marked.bam \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
-bqsr /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.GATK4.pre.recal.table


/projects/b1025/sdi0596/tools/gatk-4.1.3.0/gatk BaseRecalibrator \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
--known-sites /projects/b1025/sdi0596/whole_genome_Exome_practice/dbsnp_data/mgp.v3.snps.rsIDdbSNPv137.vcf.gz \
--use-original-qualities \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.GATK.post.recal.table

samtools index /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam


### QUALITY CONTROL OF ALIGNMENTS
java -Xmx8G -Dpicard.useLegacyParser=false \
-jar /projects/b1025/sdi0596/tools/picard.jar CollectSequencingArtifactMetrics \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I  /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam.artifacts


java -Xmx8G -Dpicard.useLegacyParser=false \
-jar /projects/b1025/sdi0596/tools/picard.jar CollectMultipleMetrics \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I  /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam.metrics


samtools idxstats /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam > /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam.idxstats


#NOW CALCULATING METRICS FOR COVERAGE CALCULATION. THIS IS TO BE DONE LATER ONCE I HAVE THE BAIT AND TARGET INTERVAL FILES
java -Xmx8G -Dpicard.useLegacyParser=false \
-jar /projects/b1025/sdi0596/tools/picard.jar CollectHsMetrics \
-SAMPLE_SIZE 100000 \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files/mm9.fa \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
-O /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam.metrics \


samtools flagstat /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam > /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.alignment.stats.txt


#### NOW calling SNV and indel
/projects/b1025/sdi0596/tools/gatk-4.1.3.0/gatk Mutect2 \
--native-pair-hmm-threads 4 \
-R /projects/b1025/sdi0596/whole_genome_Exome_practice/ref_files_GRcm38/GRCm38_68.fa \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2210079.bam \
-I /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/ERR2230866.bam \
--tumor-sample ERR2210079 --normal-sample ERR2230866 \
-O  /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/WES.m2.vcf \
-bamout /projects/b1025/sdi0596/whole_genome_Exome_practice/post_processed_file/WES.m2.bam

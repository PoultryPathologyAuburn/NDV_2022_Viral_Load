#!/bin/bash

# all the preps:
# indexing references 
# fastqc 1
# trimmomatic
# fastqc 2

#load module
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
module load trimmomatic
module load gatk/4.4.0.0
module load bwa/0.7.17
module load fastqc/0.10.1

# download ref-genome.fasta and .gtf from NCBI virus into folder ncbi-references
# .gtf for AF077761 is not available!

cd ncbi-references

# index ref-genomes:
ref_list=("ASM478661.1" "AF077761.1" "ASM283408.1")

for ref in "${ref_list[@]}"; do
	bwa index ${ref}.fasta 
	samtools faidx ${ref}.fasta
 	gatk --java-options "-Xmx1G" CreateSequenceDictionary -R ${ref}.fasta
	genome_size=`awk '{sum+=$2} END {print sum}' ${ref}.fasta.fai`
 	echo $genome_size >> genome_size.txt
done

cd ..

# fastqc raw, trimmomatic, fastqc trimmomatic:
# make sample list
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")

# make directories for the outputs
mkdir viral_fastqc1
mkdir viral_trimmomatic
mkdir viral_fastqc2

for sample in "${sample_list[@]}"; do
	echo "CURRENT FILE: $sample"
 	# check raw file with fastqc
  	fastqc -o viral_fastqc1 raw/${sample}_1.fq.gz
   	fastqc -o viral_fastqc1 raw/${sample}_2.fq.gz
    	# run trimmomatic
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar PE -phred33 raw/${sample}_1.fq.gz raw/${sample}_2.fq.gz viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
	# check paired trimmed files with fastqc
 	fastqc -o viral_fastqc2 viral_trimmomatic/f_paired_${sample}.fq.gz
   	fastqc -o viral_fastqc2 viral_trimmomatic/r_paired_${sample}.fq.gz
done

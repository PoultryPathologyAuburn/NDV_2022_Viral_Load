#!/bin/bash

# indexing references 
# fastqc
# trimmomatic

#load module
source /apps/profiles/modules_asax.sh.dyn
# module load gcc/6.2.0 bwa/0.7.12
module load samtools/1.18
# module load picard/1.79
module load gatk/4.4.0.0
module load bwa/0.7.17
module load fastqc/0.10.1

cd ncbi-references

ref_list=("ASM478661.1" "AF077761.1" "ASM283408.1")

for ref in "${ref_list[@]}"; do
	bwa index ${ref}.fasta 
	samtools faidx ${ref}.fasta
 	gatk --java-options "-Xmx1G" CreateSequenceDictionary -R ${ref}.fasta
	#java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar REFERENCE=GCF_004786615.1_ASM478661v1_genomic.fasta OUTPUT=GCF_004786615.1_ASM478661v1_genomic.dict
	genome_size=`awk '{sum+=$2} END {print sum}' ${ref}.fasta.fai`
 	echo $genome_size >> genome_size.txt
done

cd ..

# run FastQC on the raw sequences:
mkdir viral_fastqc1

for file in raw/*.fq.gz; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Perform FastQC on the file and save in directory viral_fastqc1
        fastqc -o viral_fastqc1 "$file"
    fi
done

# run trimmomatic 
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")

mkdir viral_trimmomatic

for sample in "${sample_list[@]}"; do
	echo "CURRENT FILE: $sample"
 	trimmomatic ARGS PE -phred33 raw/${sample}_1.fq.gz raw/${sample}_2.fq.gz viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz ILLUMINACLIP:raw/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        #java -jar /mnt/beegfs/apps/dmc/apps/spack_0.16.0/spack/opt/spack/linux-centos7-ivybridge/gcc-10.2.0/trimmomatic-0.39-ili2pw5eux5c4zkvobobylopjwwu7phd/bin/trimmomatic-0.39.jar PE -phred33 raw/${sample}_1.fq.gz raw/${sample}_2.fq.gz viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz ILLUMINACLIP:raw/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
done

cd ..

# run fastqc on the trimmed sequences
mkdir viral_fastqc2

for file in viral_trimmomatic/*.fq.gz; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Perform FastQC on the file and save in directory viral_fastqc2
        fastqc -o viral_fastqc2 "$file"
    fi
done

# align trimmed paired sequences against ref_genome AF077761.1, ASM478661.1, ASM283408.1  4 threads
# in case there are a lot of unpaired sequences, do the same for those
mkdir viral_PE_alignedAF
mkdir viral_PE_alignedASM47
mkdir viral_PE_alignedASM28

for sample in "${sample_list[@]}"; do
	echo "CURRENT FILE: $sample"
        bwa mem -t 4 ref_index/AF077761.1 viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz | samtools sort -o viral_PE_alignedAF/${sample}_PE.sorted.bam
	bwa mem -t 4 ref_index/ASM478661.1 viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz | samtools sort -o viral_PE_alignedASM47/${sample}_PE.sorted.bam
        bwa mem -t 4 ref_index/ASM283408.1 viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz | samtools sort -o viral_PE_alignedASM28/${sample}_PE.sorted.bam
done

mkdir viral_SE_alignedAF
mkdir viral_SE_alignedASM47
mkdir viral_SE_alignedASM28
for sample in "${sample_list[@]}"; do
	echo "CURRENT FILE: $sample"
        bwa mem -t 4 ref_index/AF077761.1 viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz | samtools sort -o viral_SE_alignedAF/${sample}_SE.sorted.bam
	bwa mem -t 4 ref_index/ASM478661.1 viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz | samtools sort -o viral_SE_alignedASM47/${sample}_SE.sorted.bam
        bwa mem -t 4 ref_index/ASM283408.1 viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz | samtools sort -o viral_SE_alignedASM28/${sample}_SE.sorted.bam
done

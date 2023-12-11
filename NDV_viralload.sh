#!/bin/bash

# indexing references 
# fastqc
# trimmomatic

#load module
source /apps/profiles/modules_dmc.sh.dyn
# module load gcc/6.2.0 bwa/0.7.12
module load samtools/1.13
# module load picard/1.79
module load gatk/4.4.0.0

mkdir ref_index
cd ref_index

ref_list=("ASM478661.1" "AF077761.1" "ASM283408.1")

for ref in "${ref_list[@]}"; do
	bwa index ../ncbi-references/${ref}.fasta 
	samtools faidx ../ncbi-references/${ref}.fasta
 	gatk --java-options "-Xmx1G" CreateSequenceDictionary -R ../ncbi-references/${ref}.fasta
	#java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar REFERENCE=GCF_004786615.1_ASM478661v1_genomic.fasta OUTPUT=GCF_004786615.1_ASM478661v1_genomic.dict
	genome_size=`awk '{sum+=$2} END {print sum}' ref_index/${ref}.fasta.fai`
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
 
sample_list=("hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")

mkdir viral_trimmomatic

for sample in "${sample_list[@]}"; do
	echo "CURRENT FILE:$sample"
        java -jar /mnt/beegfs/apps/dmc/apps/spack_0.16.0/spack/opt/spack/linux-centos7-ivybridge/gcc-10.2.0/trimmomatic-0.39-ili2pw5eux5c4zkvobobylopjwwu7phd/bin/trimmomatic-0.39.jar PE -phred33 viral_trimmomatic/${sample}_1.fq.gz viral_trimmomatic/${sample}_2.fq.gz viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/f_unpaired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz viral_trimmomatic/r_unpaired_${sample}.fq.gz ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
done

cd ..

mkdir viral_fastqc2

for file in viral_trimmomatic/*.fq.gz; do
    # Check if the file exists
    if [ -e "$file" ]; then
        # Perform FastQC on the file and save in directory viral_fastqc2
        fastqc -o viral_fastqc2 "$file"
    fi
done



#bwa mem aligns #and samtools sorts

for i in tr12v2 tr12v4 tr24v4 tr24v7 #hg12c3hg12c4hg12c5hg12v2hg12v4hg12v5hg24c2hg24c4hg24c5hg24v1hg24v4hg24v6tr12c1tr12c4tr12c5# 
          
do
        bwa mem -t 8 ndv-index output_forward_paired_${i}.fq.gz output_reverse_paired_${i}.fq.gz | samtools sort -o ${i}_output.sorted.bam

        mkdir ${i}
        mv ${i}_output.sorted.bam /${i}
done



################### Chunk 4 ################
##### GENERATE INFORMATIVE DATA ############
############################################

# Step 1: Generate flagstat output #for quality control
    #Flagstat provides statistics about the input alignment #file
    # - #1) total reads: number #of reads #in the #file
    # - 2) mapped reads: number #of reads that aligned #to the reference genome
    # - 3) Secondary (supplemental) reads: additional alignments #for a #read other than the primary alignment
    # - 4) Paired #in sequencing: number #of reads that were sequenced #in a paired-#end mode
    # - 5) Read1: number #of read1s #and Read2: number #of first #and scnd reads #in pairs
    # - 6) Properly paired: number #of #read pairs that were mapped #to the reference genome #in the proper orientation #as their mate.
    # - 7) With itself #and mate mapped: Readsthat are mapped ant their mate #is also mapped
    # - 8) Singletons: number #of reads that are mapped but their mate #is unmapped
    # - 9) With mate mapped #to a different chr: number #of reads that are mapped but their mate #is mapped #to a different chromosome
    # - 10) With mate mapped #to a different chr (mapQ>=5) number #of reads that are mapped #and their mate #is mapped #to a different chromosome, but only #if the mapping quality #is 5 #or more.

#!/usr/bin/env python3
for i in tr24v4 tr24v7 hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 
do  
    samtools flagstat ${i}_output.sorted.bam > ${i}_flagstat.txt
    chmod 777 ${i}_flagstat.txt
done


 # Step 2: Generate depth output
    #it shows a list of coordinates in the genome and number of times a read overlaps each coordinate
    #output file has three columns: chromosome (name of the reference sequence), position (position in the reference sequence where the depth is calculated), 
    #and depth (number of aligned reads that overlap the position)

for i in tr24v7 hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 
do
    samtools depth ${i}_output.sorted.bam > ${i}_depth.txt
    chmod 777 ${i}_depth.txt
done

# Step 3: Generate viral load estimate
    #output contains:
    # - no_feature: number of reads which do not overlap any feature
    # - ambiguous: number of reads which overlap more than one feature
    # - too_low_aQual: number of reads which have a mapping quality lower than the minimum quality allowed
    # - not_aligned: number of reads which do not align to the reference genome
    # - alignment_not_unique: number of reads which align to more than one location in the reference genome

for i in tr24v7 hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 
do
    python -m HTSeq.scripts.count --stranded=no -r name -f bam ${i}_output.sorted.bam genomic.gtf > ${i}_viral_load.txt
    chmod 777 ${i}_viral_load.txt
done

 # Step 4: Calculate RPM and print result

for i in tr24v7 hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 
do
    total_reads=$(awk 'NR==1{print $1}' ${i}_flagstat.txt)
    viral_load=$(awk '{sum+=$2} END {print sum}' ${i}_viral_load.txt)
    rpm=$(echo "scale=2; $viral_load / $total_reads * 1000000" | bc)
    echo "${i}: ${rpm} RPM"
done

################################################
################## CHUNKY 1 ####################
###### CREATE REF GEN INDEXES ######
################################################


#!/bin/bash
#
#load module
module load gcc/6.2.0 bwa/0.7.12
module load samtools/1.13
module load picard/1.79

bwa index GCF_004786615.1_ASM478661v1_genomic.fasta 
samtools faidx GCF_004786615.1_ASM478661v1_genomic.fasta
java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar REFERENCE=GCF_004786615.1_ASM478661v1_genomic.fasta OUTPUT=GCF_004786615.1_ASM478661v1_genomic.dict
genome_size=`awk '{sum+=$2} END {print sum}' GCF_004786615.1_ASM478661v1_genomic.fasta.fai`


################################################
################## CHUNKY 2 ####################
###### TRIM POOR SEQUENCES W TRIMMOMATICS ######
################################################

#!/bin/bash

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh

#module load trimmomatic
	
for i in hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 tr24v4 tr24v7
do 
	cd ${i}
	java -jar /mnt/beegfs/apps/dmc/apps/spack_0.16.0/spack/opt/spack/linux-centos7-ivybridge/gcc-10.2.0/trimmomatic-0.39-ili2pw5eux5c4zkvobobylopjwwu7phd/bin/trimmomatic-0.39.jar PE -phred33 ${i}_1.fq.gz ${i}_2.fq.gz output_forward_paired_${i}.fq.gz output_forward_unpaired_${i}.fq.gz output_reverse_paired_${i}.fq.gz output_reverse_unpaired_${i}.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
	cd ..
done

for i in hg12c3 hg12c4 hg12c5 hg12v2 hg12v4 hg12v5 hg24c2 hg24c4 hg24c5 hg24v1 hg24v4 hg24v6 tr12c1 tr12c4 tr12c5 tr12v2 tr12v4 tr24v4 tr24v7
do
	cd ${i}
	mv output_reverse_paired_${i}.fq.gz ..
	mv output_forward_paired_${i}.fq.gz ..
	cd ..
done

###############################################
################## CHUNKY 3 ###################
#### BWA ALIGNMENT AGAINST NDV REF GENOME #####
###############################################

#!/bin/bash

#load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2019.10
module load bwa/0.7.12
module load samtools/1.13
#
#
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

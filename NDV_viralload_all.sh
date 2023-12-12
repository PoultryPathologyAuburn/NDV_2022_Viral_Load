#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

module load gatk/4.4.0.0
module load bwa/0.7.17


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

# make sample list
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")
refdir_list=("viral_PE_alignedAF" "viral_PE_alignedASM28" "viral_PE_alignedASM47")

for refdir in "${refdir_list[@]}"; do
	for sample in "${sample_list[@]}"; do
		echo "CURRENT DIR: $refdir\nCURRENT FILE: $sample"
  		mkdir ${refdir}_results
    		results="../${refdir}_results"
		samtools flagstat ${refdir}/${sample}_output.sorted.bam > $results/${sample}_flagstat.txt
 		samtools depth ${refdir}/${sample}_output.sorted.bam > $results/${sample}_depth.txt
		python -m HTSeq.scripts.count --stranded=no -r name -f bam ${refdir}/${sample}_output.sorted.bam genomic.gtf > $results/${sample}_viral_load.txt
		total_reads=$(awk 'NR==1{print $1}' $results/${sample}_flagstat.txt)
		viral_load=$(awk '{sum+=$2} END {print sum}' $results/${sample}_viral_load.txt)
		rpm=$(echo "scale=2; $viral_load / $total_reads * 1000000" | bc)
		echo "${sample}: $total_reads totalreads $viral_load viralload ${rpm} RPM" > $results/${sample}_rpm.txt
done

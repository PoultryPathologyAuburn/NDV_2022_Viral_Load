#!/bin/bash

# stats
# counts, total reads, viral load, RPM

source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18
module load htseq/2.4.0

# make sample list
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")
ref_list=("AF077761.1" "ASM283408.1" "ASM478661.1")

for ref in "${ref_list[@]}"; do
	mkdir ${ref}_results
	for sample in "${sample_list[@]}"; do
		echo "CURRENT REF: $ref\nCURRENT FILE: $sample"
		samtools flagstat ${ref}/${sample}_sorted.bam > ${ref}_results/${sample}_flagstat.txt
 		samtools depth ${ref}/${sample}_sorted.bam > ${ref}_results/${sample}_depth.txt
		python -m HTSeq.scripts.count --stranded=no -r name -f bam ${ref}/${sample}_sorted.bam ncbi-references/${ref}_genomic.gtf > ${ref}_results/${sample}_counts.txt
		total_reads=$(awk 'NR==1{print $1}' ${ref}_results/${sample}_flagstat.txt)
		viral_load=$(awk '{sum+=$2} END {print sum}' ${ref}_results/${sample}_counts.txt)
		rpm=$(echo "scale=2; $viral_load / $total_reads * 1000000" | bc)
		echo "${sample}: $total_reads totalreads $viral_load viralload ${rpm} RPM" >> ${ref}_results/${ref}_stats.txt
	done
done

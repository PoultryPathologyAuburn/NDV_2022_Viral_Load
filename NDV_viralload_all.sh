#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

# stats

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
		python -m HTSeq.scripts.count --stranded=no -r name -f bam ${refdir}/${sample}_output.sorted.bam ncbi-references/${refdir}_genomic.gtf > $results/${sample}_viral_load.txt
		total_reads=$(awk 'NR==1{print $1}' $results/${sample}_flagstat.txt)
		viral_load=$(awk '{sum+=$2} END {print sum}' $results/${sample}_viral_load.txt)
		rpm=$(echo "scale=2; $viral_load / $total_reads * 1000000" | bc)
		echo "${sample}: $total_reads totalreads $viral_load viralload ${rpm} RPM" >> stats.txt
	done
done

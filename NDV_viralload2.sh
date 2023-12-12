#!/bin/bash

# make sample list
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")
ref_list=("ASM478661.1" "AF077761.1" "ASM283408.1")

for ref in "${ref_list[@]}"; do
  mkdir ${ref}
  for sample in "${sample_list[@]}"; do
    echo "CURRENT ITEM $ref/$sample"
    bwa mem -t 4 ncbi-references/${ref}.fasta viral_trimmomatic/${sample}_1.fq.gz viral_trimmomatic/${sample}_2.fq.gz > ${ref}/${sample}.sam
    samtools view -S -b ${ref}/${sample}.sam | samtools sort -o ${ref}/${sample}_sorted.bam - && samtools index ${ref}/${sample}_sorted.bam
    echo "${ref}/${sample} $(samtools view -c -F 4 ${ref}/${sample}_sorted.bam)" >> ${ref}/${ref}_mapped_counts.txt
    htseq-count -f bam -r pos -s no -i gene_id ${ref}/${sample}_sorted.bam ncbi-references/${ref}.gtf
  done
done   

#samtools view -S -b LaSota.sam > LaSota.bam
#samtools sort LaSota.bam -o LaSota_sorted.bam
#samtools index LaSota_sorted.bam
#samtools view -c -F 4 LaSota_sorted.bam
#htseq-count -f bam -r pos -s no -i gene_id sorted_aligned_reads.bam reference_annotation.gtf


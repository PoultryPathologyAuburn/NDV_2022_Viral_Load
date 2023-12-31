#!/bin/bash

# mapping reads against references
# convert, sort, index
# count hits

# load modules
source /apps/profiles/modules_asax.sh.dyn
module load bwa/0.7.12
module load samtools/1.18
module load htseq/2.0.4

# make sample list
sample_list=("LaSota" "hg12c3" "hg12c4" "hg12c5" "hg12v2" "hg12v4" "hg12v5" "hg24c2" "hg24c4" "hg24c5" "hg24v1" "hg24v4" "hg24v6" "tr12c1" "tr12c4" "tr12c5" "tr12v2" "tr12v4" "tr24v4" "tr24v7")
ref_list=("ASM478661.1" "AF077761.1" "ASM283408.1")

for ref in "${ref_list[@]}"; do
  mkdir ${ref}
  for sample in "${sample_list[@]}"; do
    echo "CURRENT ITEM $ref/$sample"
    # map
    bwa mem -t 4 ncbi-references/${ref}.fasta viral_trimmomatic/f_paired_${sample}.fq.gz viral_trimmomatic/r_paired_${sample}.fq.gz > ${ref}/${sample}.sam
    # convert sam to bam, sort, index
    samtools view -S -b ${ref}/${sample}.sam | samtools sort -o ${ref}/${sample}_sorted.bam - && samtools index ${ref}/${sample}_sorted.bam
    # samtools count
    samtools_counts=$(samtools view -c -F 4 ${ref}/${sample}_sorted.bam)
    # in case "htseq_counts" does not work, sub it with "python -m HTSeq.scripts.count"
    # htseq_counts
    htseq_counts=$(htseq-count -f bam -r pos -s no -i gene_id ${ref}/${sample}_sorted.bam ncbi-references/${ref}.gtf)
    echo "${ref}/${sample} $samtools_counts $htseq_counts" >> ${ref}/${ref}_counts.txt
  done
done   


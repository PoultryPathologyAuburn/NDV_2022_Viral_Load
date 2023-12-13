# NDV_2022_Viral_Load
script for indexing, fastqc, trimming: <a href="NDV_prep.sh">NDV_prep.sh</a> <br>
script for mapping against reference genome: <a href="NDV_map.sh">NDV_map.sh</a> <br>
script for stats: <a href="NDV_stat.sh">NDV_stat.sh</a>
<details>
  <ul>
    <li>SPF chickens, inoculated with NDV, samples taken after 12, 24, 48 hours, from Harderian gland, trachea, spleen, cecal tonsils</li>
    <li>NCBI virus reference genomes: AF077761.1.fna (LaSota 1999), ASM283408.1.fna (Chinese ducks 2018), ASM478661.1.fna (Ireland chickens 2023)</li>
    <li>positive control: LaSota sequences from Kyriakis lab</li>
  </ul>
</details>

## 1. Indexing reference genomes: 
<ul>
    <li>bwa index</li> 
    <li>samtools faidx</li>
    <li>GATK CreateSequenceDictionary</li>
    <li>NCBI reference genomes: AF077761.1.fna (LaSota 1999), ASM283408.1.fna (Chinese ducks 2018), ASM478661.1.fna (Ireland chickens 2023)</li></ul>

## 2. Quality check with FastQC
<ul>
    <li>output: viral_fastqc1</li> </ul>
  
## 3. Trimming with Trimmomatic
<ul>
    <li>output: viral_trimmomatic</li> </ul>

## 4. Quality check with FastQC
<ul>
    <li>output: viral_fastqc2 (not done!)</li> </ul>
  
## 5. Alignment against all three references and counting hits: 
<ul>
    <li>bwa mem</li>
    <li>samtools sort</li>
    <li>samtools view</li>
    <li>htseq-count</li>
  <li>for all three references, for trimmed sequences, for raw sequences</li>
    <li>output _sorted.bam: for trimmed the output folders have the reference-name, for raw reference_raw</li></ul>
    
## 6. Statistics: 
<ul>
    <li>samtools flagstat</li>
    <li>samtools depth</li>
    <li>HTSeq (viral load)</li>
    <li>rpm: viral load/total reads</li>
    <li>output: references_results</li></ul>

## group permission for shared folders:
<kbd>chmod -R +777 foldername</kbd>

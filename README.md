# NDV_2022_Viral_Load
script: <a href="NDV_viralload.sh">NDV_viralload.sh</a>
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
    <li>Picard CreateSequenceDictionary</li>
    <li>NCBI reference genomes: AF077761.1.fna (LaSota 1999), ASM283408.1.fna (Chinese ducks 2018), ASM478661.1.fna (Ireland chickens 2023)</li></ul>

## 2. Quality check with FastQC
  
## 3. Trimming with Trimmomatic

## 4. Quality check with FastQC
  
## 5. Alignment against reference: 
<ul>
    <li>bwa mem</li>
    <li>samtools sort</li></ul>
    
## 6. Statistics: 
<ul>
    <li>samtools flagstat</li>
    <li>samtools depth</li>
    <li>HTSeq (viral load)</li>
    <li>viral load/total reads</li></ul>
